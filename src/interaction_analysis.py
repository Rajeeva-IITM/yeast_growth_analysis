"""Partial-dependence and interaction analysis for the bloom_256 encoding runs.
"""

import re
from collections import Counter
from pathlib import Path
from typing import Optional

import hydra
import lightgbm as lgb
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

matplotlib.rcParams["svg.fonttype"] = "none"  # keep panel text editable in the per-panel SVGs
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import polars as pl  # noqa: E402
import polars.selectors as cs  # noqa: E402
from dotenv import load_dotenv  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from numpy import ndarray  # noqa: E402
from omegaconf import DictConfig  # noqa: E402
from PyALE._src.ALE_1D import aleplot_1D_discrete  # noqa: E402
from PyALE._src.ALE_2D import aleplot_2D_continuous  # noqa: E402
from rich.console import Console  # noqa: E402
from sklearn.base import BaseEstimator, RegressorMixin  # noqa: E402
from sklearn.inspection import partial_dependence  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import train_test_split  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402
from utils import get_model, get_model_paths  # noqa: E402

load_dotenv()

console = Console()

META_COLUMNS = ("Strain", "Condition", "Phenotype")


# --------------------------------------------------------------------------------------
# Data / model plumbing
# --------------------------------------------------------------------------------------


def seed_from_run_dir(model_path: Path) -> int:
    """Recover the training seed from a run directory named ``<Prefix>_<Name>_<seed>_<suffix>``.

    The seed is needed to reproduce the exact ``train_test_split`` the model was fitted
    under, so that the StandardScaler used here matches the one the model saw.

    Parameters
    ----------
    model_path : Path
        Path to ``Boosting.pkl`` inside the run directory.

    Returns
    -------
    int
        The seed encoded in the directory name.
    """
    parts = model_path.parent.name.split("_")
    for token in reversed(parts):
        if token.isdigit():
            return int(token)
    raise ValueError(f"no seed found in run directory name {model_path.parent.name!r}")


def load_dataset(data_path: str) -> tuple[ndarray, ndarray, ndarray, list[str]]:
    """Load one encoding dataset in the exact column order ``get_data`` produces.

    ``utils.get_data(run_type="full")`` selects ``cs.starts_with("Y")`` then
    ``cs.contains("latent")``; the model's feature indices follow that order, so the
    feature names must be built the same way.

    Parameters
    ----------
    data_path : str
        Path to the ``.feather`` file.

    Returns
    -------
    tuple[ndarray, ndarray, ndarray, list[str]]
        Features, phenotype, condition labels, and feature names.
    """
    df = pl.read_ipc(data_path)
    features = df.select(cs.starts_with("Y"), cs.contains("latent"))
    names = features.columns
    X = features.to_numpy()
    y = df["Phenotype"].to_numpy()
    conditions = df["Condition"].to_numpy()
    del df, features
    return X, y, conditions, names


def fit_run_scaler(
    X: ndarray, y: ndarray, seed: int, test_frac: float, verify: bool, model: lgb.Booster
) -> tuple[StandardScaler, Optional[float]]:
    """Reproduce a run's train/test split and fit its StandardScaler on the train half.

    ``tune_model.py`` splits with ``train_test_split(test_size=test_frac,
    random_state=seed)`` and fits the scaler on the training portion only. Splitting an
    index array reproduces the identical partition, which avoids materialising a second
    copy of a multi-GB feature matrix.

    Parameters
    ----------
    X, y : ndarray
        Full feature matrix and target.
    seed : int
        The run's seed.
    test_frac : float
        Test fraction, from ``conf.testing.test_frac``.
    verify : bool
        If True, also score the model on the held-out half and return its AUC.
    model : lgb.Booster
        Used only for the verification score.

    Returns
    -------
    tuple[StandardScaler, Optional[float]]
        The fitted scaler and, when ``verify``, the held-out ROC AUC.
    """
    index = np.arange(X.shape[0])
    train_idx, test_idx = train_test_split(
        index, test_size=test_frac, random_state=seed
    )

    scaler = StandardScaler().fit(X[train_idx])

    auc = None
    if verify:
        preds = model.predict(scaler.transform(X[test_idx]))
        auc = float(roc_auc_score(y[test_idx], preds))

    return scaler, auc


def stratified_sample(
    conditions: ndarray, n_sample: int, rng: np.random.Generator
) -> ndarray:
    """Pick row indices spread evenly over conditions.

    A uniform random sample would leave rare conditions with too few rows to estimate a
    condition-specific partial-dependence curve. Sampling per condition keeps every
    condition represented.

    Parameters
    ----------
    conditions : ndarray
        Per-row condition labels.
    n_sample : int
        Target number of rows.
    rng : np.random.Generator
        Source of randomness.

    Returns
    -------
    ndarray
        Sorted row indices.
    """
    unique = np.unique(conditions)
    per_condition = max(1, n_sample // len(unique))

    picked = []
    for condition in unique:
        rows = np.flatnonzero(conditions == condition)
        take = min(per_condition, len(rows))
        picked.append(rng.choice(rows, size=take, replace=False))

    return np.sort(np.concatenate(picked))


# --------------------------------------------------------------------------------------
# Partial dependence
# --------------------------------------------------------------------------------------


def build_grid(values: ndarray, max_grid: int) -> ndarray:
    """Grid points for one feature, taken from values that actually occur.

    Called on the background sample, so the grid matches the one
    ``sklearn.inspection.partial_dependence`` derives from the same rows.

    Gene columns hold {0, 1, 2, 3} (max/sum) or small counts, so their full support is a
    handful of points. Latent columns hold one value per condition -- at most 39 -- so they
    are also enumerable. Only if a feature exceeds ``max_grid`` distinct values does this
    fall back to quantiles, which keeps every grid point on-manifold.

    Parameters
    ----------
    values : ndarray
        Observed values of the feature.
    max_grid : int
        Maximum number of grid points.

    Returns
    -------
    ndarray
        Ascending grid values, in the feature's original units.
    """
    unique = np.unique(values)
    if len(unique) <= max_grid:
        return unique
    return np.unique(np.quantile(values, np.linspace(0, 1, max_grid)))


class BoosterEstimator(BaseEstimator, RegressorMixin):
    """Minimal sklearn estimator wrapping a raw ``lgb.Booster``.

    ``sklearn.inspection.partial_dependence`` needs an estimator, and a booster unpickled
    from a run directory is not one. Registering as a regressor keeps sklearn's
    ``response_method`` on ``predict``, so partial dependence is computed in raw-score
    space, where contributions are additive -- the same space the hand-rolled sweep used.

    Parameters
    ----------
    booster : lgb.Booster
        The fitted model to wrap.
    """

    def __init__(self, booster: Optional[lgb.Booster] = None):
        self.booster = booster

    def fit(self, X: ndarray, y: Optional[ndarray] = None) -> "BoosterEstimator":
        """Record the input width sklearn validates against. The booster must be already fit."""
        self.n_features_in_ = X.shape[1]
        self.is_fitted_ = True
        return self

    def predict(self, X: ndarray) -> ndarray:
        """Raw scores, matching what the rest of this module works in.

        Coerced with ``np.asarray`` because PyALE passes a named-column DataFrame while the
        booster was trained on a bare array; LightGBM would otherwise reject the feature
        names.
        """
        return self.booster.predict(np.asarray(X), raw_score=True)


def _check_grid(returned: ndarray, expected: ndarray, label: str) -> None:
    """Fail loudly if sklearn did not use the grid we asked for.

    sklearn picks its own grid: when a feature has more distinct values than
    ``grid_resolution`` it falls back to ``np.linspace(min, max, n)``, which places points
    at values that never occur in the data. For the latent block -- only 18-39 distinct
    chemistry vectors exist -- that would be pure extrapolation, exactly the artefact this
    module is built to avoid. ``build_grid`` is the specification; this is the tripwire.

    Parameters
    ----------
    returned : ndarray
        The grid sklearn actually used.
    expected : ndarray
        The grid ``build_grid`` asked for.
    label : str
        Feature identifier, for the error message.

    Raises
    ------
    ValueError
        If the grids differ.
    """
    same_length = len(returned) == len(expected)
    if same_length and np.allclose(np.sort(returned), np.sort(expected)):
        return

    detail = (
        f"grid sizes differ ({len(returned)} vs {len(expected)} requested)"
        if not same_length
        else f"both grids have {len(returned)} points but the values differ"
    )
    raise ValueError(
        f"sklearn did not use the requested grid for {label}: {detail}. The feature has "
        "more distinct values than analysis.max_grid, so sklearn fell back to a linspace "
        "grid whose points do not occur in the data -- partial dependence there would be "
        "extrapolation. Raise analysis.max_grid above the feature's cardinality (37 is "
        "the highest seen in these datasets), or set analysis.pdp_backend=builtin, which "
        "grids on quantiles and stays on-manifold."
    )


def _pdp_sweep_builtin(
    model: lgb.Booster, X_scaled: ndarray, column: int, grid_scaled: ndarray
) -> ndarray:
    """Hand-rolled equivalent of :func:`pdp_sweep`, kept for cross-checking.

    Verified bit-identical to the sklearn path (``max|diff| = 0.0``). Reachable via
    ``analysis.pdp_backend=builtin`` so the equivalence can be re-demonstrated rather than
    taken on trust.
    """
    original = X_scaled[:, column].copy()
    surface = np.empty((len(grid_scaled), X_scaled.shape[0]))

    for i, value in enumerate(grid_scaled):
        X_scaled[:, column] = value
        surface[i] = model.predict(X_scaled, raw_score=True)

    X_scaled[:, column] = original
    return surface


def _pdp_2d_builtin(
    model: lgb.Booster,
    X_scaled: ndarray,
    col_a: int,
    col_b: int,
    grid_a: ndarray,
    grid_b: ndarray,
) -> ndarray:
    """Hand-rolled equivalent of :func:`pdp_2d`, kept for cross-checking."""
    original_a = X_scaled[:, col_a].copy()
    original_b = X_scaled[:, col_b].copy()
    out = np.empty((len(grid_a), len(grid_b)))

    for i, value_a in enumerate(grid_a):
        X_scaled[:, col_a] = value_a
        for j, value_b in enumerate(grid_b):
            X_scaled[:, col_b] = value_b
            out[i, j] = model.predict(X_scaled, raw_score=True).mean()

    X_scaled[:, col_a] = original_a
    X_scaled[:, col_b] = original_b
    return out


def pdp_sweep(
    model: lgb.Booster,
    X_scaled: ndarray,
    column: int,
    grid_scaled: ndarray,
    max_grid: int = 40,
    backend: str = "sklearn",
) -> ndarray:
    """Predict over a feature's grid, holding every other feature at its observed value.

    Returns the full per-row prediction surface rather than the column means, because the
    caller needs it three ways: averaged (marginal PDP), grouped by condition (G x E) and
    per row (ICE). Computing it once and reducing it three ways is what makes the G x E
    analysis free.

    Parameters
    ----------
    model : lgb.Booster
        The model to probe.
    X_scaled : ndarray
        Background rows, already in the model's scaled input space. Modified in place and
        restored before returning.
    column : int
        Index of the feature to sweep.
    grid_scaled : ndarray
        Grid values, in scaled units, as chosen by :func:`build_grid`.
    max_grid : int
        Passed to sklearn as ``grid_resolution``; must exceed the feature's cardinality or
        :func:`_check_grid` raises.
    backend : str
        ``"sklearn"`` or ``"builtin"``.

    Returns
    -------
    ndarray
        Shape ``(n_grid, n_rows)`` of raw-score predictions.
    """
    if backend == "builtin":
        return _pdp_sweep_builtin(model, X_scaled, column, grid_scaled)

    estimator = BoosterEstimator(model).fit(X_scaled)
    result = partial_dependence(
        estimator,
        X_scaled,
        [column],
        method="brute",
        percentiles=(0, 1),
        grid_resolution=max_grid,
        kind="both",
    )
    _check_grid(result["grid_values"][0], grid_scaled, f"column {column}")

    # kind="both" hands back the per-row surface as well as its mean, which is what keeps
    # the G x E analysis free: one sweep, reduced three ways downstream.
    return np.asarray(result["individual"])[0].T


def pdp_2d(
    model: lgb.Booster,
    X_scaled: ndarray,
    col_a: int,
    col_b: int,
    grid_a: ndarray,
    grid_b: ndarray,
    max_grid: int = 40,
    backend: str = "sklearn",
) -> ndarray:
    """Two-dimensional partial dependence over a pair of features.

    Parameters
    ----------
    model : lgb.Booster
        The model to probe.
    X_scaled : ndarray
        Background rows in scaled space. Modified in place and restored.
    col_a, col_b : int
        Feature indices.
    grid_a, grid_b : ndarray
        Grids in scaled units, as chosen by :func:`build_grid`.
    max_grid : int
        Passed to sklearn as ``grid_resolution``.
    backend : str
        ``"sklearn"`` or ``"builtin"``.

    Returns
    -------
    ndarray
        Shape ``(len(grid_a), len(grid_b))`` of mean raw scores.
    """
    if backend == "builtin":
        return _pdp_2d_builtin(model, X_scaled, col_a, col_b, grid_a, grid_b)

    estimator = BoosterEstimator(model).fit(X_scaled)
    result = partial_dependence(
        estimator,
        X_scaled,
        [(col_a, col_b)],
        method="brute",
        percentiles=(0, 1),
        grid_resolution=max_grid,
        kind="average",
    )
    _check_grid(result["grid_values"][0], grid_a, f"column {col_a}")
    _check_grid(result["grid_values"][1], grid_b, f"column {col_b}")

    return np.asarray(result["average"])[0]


def friedman_h2(pd_ab: ndarray, pd_a: ndarray, pd_b: ndarray) -> float:
    """Friedman's H^2 for one feature pair.

    Centres the joint and marginal partial-dependence surfaces, then reports the share of
    the joint surface's variance that the two marginals fail to explain. 0 means the pair
    is perfectly additive; 1 means the joint effect is pure interaction.

    Parameters
    ----------
    pd_ab : ndarray
        Joint partial dependence, shape ``(n_a, n_b)``.
    pd_a, pd_b : ndarray
        Marginal partial dependence for each feature.

    Returns
    -------
    float
        H^2 in [0, 1], or ``nan`` when the joint surface is flat.
    """
    joint = pd_ab - pd_ab.mean()
    marginal_a = (pd_a - pd_a.mean())[:, None]
    marginal_b = (pd_b - pd_b.mean())[None, :]

    residual = joint - marginal_a - marginal_b
    denominator = float((joint**2).sum())

    if denominator <= 0:
        return float("nan")
    return float((residual**2).sum() / denominator)


def gxe_interaction(surface: ndarray, condition_codes: ndarray, n_conditions: int) -> float:
    """Two-way interaction ratio between a gene's dosage and the condition.

    Builds the ``(n_grid, n_conditions)`` table of condition-specific partial dependence
    and reports the share of its variance left after removing the gene main effect and the
    condition main effect. This is Friedman's H^2 with the condition treated as a single
    categorical factor, which is the correct way to handle chemistry here: the 256 latent
    dimensions are a deterministic function of the condition, so varying them
    independently would be off-manifold.

    Parameters
    ----------
    surface : ndarray
        Prediction surface from :func:`pdp_sweep`, shape ``(n_grid, n_rows)``.
    condition_codes : ndarray
        Integer condition code per background row.
    n_conditions : int
        Number of distinct conditions.

    Returns
    -------
    float
        Interaction ratio in [0, 1], or ``nan`` when the gene has no effect at all.
    """
    table = np.empty((surface.shape[0], n_conditions))
    for code in range(n_conditions):
        table[:, code] = surface[:, condition_codes == code].mean(axis=1)

    grand = table.mean()
    gene_effect = table.mean(axis=1, keepdims=True) - grand
    condition_effect = table.mean(axis=0, keepdims=True) - grand

    centred = table - grand
    residual = centred - gene_effect - condition_effect
    denominator = float((centred**2).sum())

    if denominator <= 0:
        return float("nan")
    return float((residual**2).sum() / denominator)



# --------------------------------------------------------------------------------------
# ALE (PyALE)
# --------------------------------------------------------------------------------------


def weighted_var(values: ndarray, weights: ndarray) -> float:
    """Variance of an effect surface under the observed distribution of its features.

    ``Var(f_j(X_j))`` in the H^2_ALE formula is a variance over the *data* distribution,
    not over grid points. Weighting by observed cell counts is also what gives this
    statistic its resistance to linkage disequilibrium: a genotype combination that never
    occurs carries zero weight, so the model's arbitrary response there cannot enter the
    number. ALE's own binning cannot do this for binary features -- one interval means one
    cell, and that cell always contains every row.

    Parameters
    ----------
    values : ndarray
        Effect values, flattened.
    weights : ndarray
        Observed counts aligned with ``values``.

    Returns
    -------
    float
        Frequency-weighted variance, or 0.0 if nothing is observed.
    """
    total = weights.sum()
    if total <= 0:
        return 0.0
    w = weights / total
    mean = float((w * values).sum())
    return float((w * (values - mean) ** 2).sum())


def ale_1d(frame: pd.DataFrame, estimator: BoosterEstimator, name: str) -> pd.DataFrame:
    """First-order ALE for one gene column, via PyALE's discrete path.

    The discrete variant is correct here: gene columns hold a handful of integer levels,
    and PyALE's continuous path would quantile-bin them pointlessly.

    Returns
    -------
    pd.DataFrame
        Indexed by grid value, with ``eff`` (centred effect) and ``size`` (bin count).
    """
    return aleplot_1D_discrete(X=frame, model=estimator, feature=name, include_CI=False)


def ale_2d(
    frame: pd.DataFrame,
    estimator: BoosterEstimator,
    name_a: str,
    name_b: str,
    grid_size: int,
) -> pd.DataFrame:
    """Second-order ALE surface for a pair of genes.

    ``impute_empty_cells=False`` matters: PyALE's default fills unobserved cells by kd-tree
    nearest neighbour, which would invent effects for exactly the genotype combinations
    linkage disequilibrium makes impossible.

    Returns
    -------
    pd.DataFrame
        Index = grid of ``name_a``, columns = grid of ``name_b``, values = accumulated
        centred second-order effect.
    """
    return aleplot_2D_continuous(
        X=frame,
        model=estimator,
        features=[name_a, name_b],
        grid_size=grid_size,
        impute_empty_cells=False,
    )


def joint_weights(
    raw_a: ndarray, raw_b: ndarray, grid_a: ndarray, grid_b: ndarray
) -> ndarray:
    """Observed row counts for each cell of an ALE surface.

    Rows are assigned to the nearest grid coordinate in each dimension, so the weights line
    up with the surface PyALE returns.

    Returns
    -------
    ndarray
        Shape ``(len(grid_a), len(grid_b))`` of counts.
    """
    ia = np.abs(raw_a[:, None] - grid_a[None, :]).argmin(axis=1)
    ib = np.abs(raw_b[:, None] - grid_b[None, :]).argmin(axis=1)
    counts = np.zeros((len(grid_a), len(grid_b)))
    np.add.at(counts, (ia, ib), 1.0)
    return counts


def h2_ale(
    surface: pd.DataFrame, eff_a: pd.DataFrame, eff_b: pd.DataFrame, weights: ndarray
) -> float:
    """Normalised ALE interaction ratio.

    ``H2 = Var(f_jk) / (Var(f_j) + Var(f_k) + Var(f_jk))``, every variance taken under the
    observed distribution. Note this is **not** on the same scale as Friedman's H^2, which
    divides by the joint surface variance alone -- this one is systematically smaller
    because the denominator carries both main effects. Both are reported so they are not
    mistaken for one another.

    Parameters
    ----------
    surface : pd.DataFrame
        Second-order ALE surface from :func:`ale_2d`.
    eff_a, eff_b : pd.DataFrame
        First-order ALE tables from :func:`ale_1d`, carrying ``eff`` and ``size``.
    weights : ndarray
        Observed joint cell counts from :func:`joint_weights`.

    Returns
    -------
    float
        Ratio in [0, 1], or ``nan`` if every effect is flat.
    """
    values = np.nan_to_num(surface.to_numpy(), nan=0.0)
    var_joint = weighted_var(values.ravel(), weights.ravel())
    var_a = weighted_var(eff_a["eff"].to_numpy(), eff_a["size"].to_numpy().astype(float))
    var_b = weighted_var(eff_b["eff"].to_numpy(), eff_b["size"].to_numpy().astype(float))

    denominator = var_a + var_b + var_joint
    if denominator <= 0:
        return float("nan")
    return float(var_joint / denominator)

# --------------------------------------------------------------------------------------
# Feature ranking
# --------------------------------------------------------------------------------------


def permutation_importance(
    model: lgb.Booster,
    X_scaled: ndarray,
    y: ndarray,
    columns: ndarray,
    n_repeats: int,
    rng: np.random.Generator,
) -> ndarray:
    """AUC drop when each candidate feature is shuffled.

    Split gain is not trusted as the final ranking here: it records where trees split and
    how much training loss that reduced, which describes how the model was *fitted* rather
    than how it behaves. Permutation importance interrogates the fitted model the same way
    partial dependence does -- perturb one feature, measure the change in prediction -- so
    the ranking and the curves it selects rest on the same assumption.

    Parameters
    ----------
    model : lgb.Booster
        The model to probe.
    X_scaled : ndarray
        Background rows in scaled space. Restored after each feature.
    y : ndarray
        Labels for the background rows.
    columns : ndarray
        Candidate feature indices.
    n_repeats : int
        Shuffles per feature.
    rng : np.random.Generator
        Source of randomness.

    Returns
    -------
    ndarray
        Mean AUC drop per candidate, aligned with ``columns``.
    """
    baseline = roc_auc_score(y, model.predict(X_scaled))
    drops = np.zeros(len(columns))

    for i, column in enumerate(columns):
        original = X_scaled[:, column].copy()
        total = 0.0
        for _ in range(n_repeats):
            X_scaled[:, column] = rng.permutation(original)
            total += baseline - roc_auc_score(y, model.predict(X_scaled))
        X_scaled[:, column] = original
        drops[i] = total / n_repeats

    return drops


def is_gene(name: str) -> bool:
    """True for genotype columns, False for the latent chemistry block."""
    return not name.startswith("latent")


def load_shap_ranking(
    shap_path: Optional[str],
    name: str,
    model_type: str,
    feature_names: list[str],
) -> tuple[Optional[ndarray], Optional[ndarray]]:
    """Seed-averaged mean |SHAP| from ``src/shap_analysis.py``, aligned to this dataset.

    Used only to widen the candidate pool that permutation importance re-scores. No value
    returned here enters a reported statistic. Split gain and mean |SHAP| disagree about
    what matters -- one describes how the model was fitted, the other attributes its
    output observationally -- and each nominates 31-43 features per run that the other
    misses, of which 1-8 survive into the reported top 50. Taking the union stops either
    criterion deciding on its own what permutation importance is ever allowed to measure.

    The csv is written unsorted and in shap's own column order, which is not the order
    :func:`load_dataset` produces, so the join is by feature name.

    Parameters
    ----------
    shap_path : str or None
        Directory holding ``shap_{name}_{model_type}.csv``. A missing directory or file is
        not an error: ``sum`` / ``count`` / ``sigma_1.0`` have no SHAP output yet and fall
        back to the gain-only pool.
    name : str
        Dataset name, e.g. ``Bloom2013``.
    model_type : str
        Model type used in the filename, e.g. ``Boosting``.
    feature_names : list[str]
        Column names in this dataset's order.

    Returns
    -------
    tuple[ndarray or None, ndarray or None]
        Column indices sorted by descending mean |SHAP|, and the values themselves
        aligned to ``feature_names``. ``(None, None)`` when no file was found.

    Raises
    ------
    ValueError
        If the csv's feature set differs from ``feature_names``. That means ``shap_path``
        points at a different encoding or threshold, and quietly analysing a partial join
        would be worse than stopping.
    """
    if shap_path is None:
        return None, None

    path = Path(shap_path) / f"shap_{name}_{model_type}.csv"
    if not path.is_file():
        console.log(
            f"[{name}] [yellow]no SHAP file at {path}[/yellow]; "
            "candidate pool falls back to gain only"
        )
        return None, None

    table = pl.read_csv(path)
    lookup = dict(zip(table["Feature"].to_list(), table["Value"].to_list()))

    missing = [n for n in feature_names if n not in lookup]
    extra = len(lookup) - (len(feature_names) - len(missing))
    if missing or extra:
        raise ValueError(
            f"{path} does not describe this dataset: {len(missing)} of "
            f"{len(feature_names)} columns absent from the csv and {extra} csv rows "
            f"unmatched (e.g. {missing[:3]}). Check that shap_path, encoding and alpha "
            "agree with the data being analysed."
        )

    values = np.array([lookup[n] for n in feature_names], dtype=float)
    return np.argsort(-values), values


# --------------------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------------------


def plot_pdp_page(
    axis: plt.Axes,
    grid: ndarray,
    curves: ndarray,
    ice: Optional[ndarray],
    name: str,
    index: Optional[int] = None,
) -> None:
    """Draw one PDP panel: per-seed curves, their mean, and ICE lines behind.

    Discrete features (genes take at most a handful of values) are drawn as markers with a
    spread bar rather than a line, because a continuous curve through four points implies
    a smoothness the model does not have.

    Parameters
    ----------
    axis : plt.Axes
        Target axes.
    grid : ndarray
        Grid values in original units.
    curves : ndarray
        Shape ``(n_seeds, n_grid)`` of per-seed partial dependence.
    ice : Optional[ndarray]
        Shape ``(n_ice, n_grid)`` of individual curves from the first seed.
    name : str
        Feature name, used as the panel title.
    """
    discrete = len(grid) <= 6

    if ice is not None and len(ice):
        centred = ice - ice.mean(axis=1, keepdims=True)
        offset = curves.mean(axis=0).mean()
        for row in centred:
            axis.plot(grid, row + offset, color="0.75", linewidth=0.4, alpha=0.5, zorder=1)

    mean = curves.mean(axis=0)
    spread = curves.std(axis=0)

    if discrete:
        axis.errorbar(
            grid, mean, yerr=spread, fmt="o-", color="#1f77b4",
            capsize=3, linewidth=1.5, markersize=5, zorder=3,
        )
    else:
        axis.fill_between(
            grid, mean - spread, mean + spread, color="#1f77b4", alpha=0.25, zorder=2
        )
        axis.plot(grid, mean, color="#1f77b4", linewidth=1.5, zorder=3)

    axis.set_title(name if index is None else f"({index}) {name}", fontsize=9)
    axis.set_xlabel("feature value", fontsize=7)
    axis.set_ylabel("partial dependence (raw score)", fontsize=7)
    axis.tick_params(labelsize=6)
    axis.grid(alpha=0.2)


def write_pdp_pdf(
    path: Path,
    names: list[str],
    grids: dict[str, ndarray],
    curves: dict[str, ndarray],
    ice: dict[str, ndarray],
    title: str,
) -> None:
    """Write the top-N partial-dependence panels, six to a page."""
    per_page = 6
    pages = max(1, -(-len(names) // per_page))
    with PdfPages(path) as pdf:
        for start in range(0, len(names), per_page):
            chunk = names[start : start + per_page]
            figure, axes = plt.subplots(2, 3, figsize=(13, 7.5))
            for offset, (axis, name) in enumerate(zip(axes.ravel(), chunk)):
                plot_pdp_page(
                    axis, grids[name], curves[name], ice.get(name), name,
                    index=start + offset + 1,
                )
            for axis in axes.ravel()[len(chunk) :]:
                axis.axis("off")
            figure.suptitle(
                f"{title}  —  panels {start + 1}-{start + len(chunk)} of {len(names)} "
                f"(page {start // per_page + 1} of {pages})",
                fontsize=11,
            )
            figure.tight_layout(rect=(0, 0, 1, 0.96))
            pdf.savefig(figure)
            plt.close(figure)


def plot_gxe_panel(
    axis,
    gene: str,
    grid: ndarray,
    table: ndarray,
    spread: ndarray,
    score: float,
    condition_names: list[str],
    colours: ndarray,
    index: Optional[int] = None,
) -> None:
    """Draw one gene's genotype x environment panel onto an existing axis.

    Shared by :func:`write_gxe_pdf` and :func:`write_gxe_panels`, so the multi-page PDF and
    the per-gene files cannot drift apart.
    """
    centred = table - table.mean(axis=0, keepdims=True)
    for code, name in enumerate(condition_names):
        axis.fill_between(
            grid,
            centred[:, code] - spread[:, code],
            centred[:, code] + spread[:, code],
            color=colours[code], alpha=0.12, linewidth=0,
        )
        axis.plot(
            grid, centred[:, code], color=colours[code],
            linewidth=1.0, marker="o", markersize=3, label=name,
        )
    axis.axhline(0.0, color="0.4", linewidth=0.7, linestyle="--")
    label = gene if index is None else f"({index}) {gene}"
    axis.set_title(f"{label}   H²(gene × condition) = {score:.3f}", fontsize=9)
    axis.set_xlabel("mutational-status", fontsize=7)
    axis.set_ylabel("condition-centred partial dependence", fontsize=7)
    axis.tick_params(labelsize=6)
    axis.grid(alpha=0.2)


def write_gxe_pdf(
    path: Path,
    genes: list[str],
    grids: dict[str, ndarray],
    tables: dict[str, ndarray],
    spreads: dict[str, ndarray],
    scores: dict[str, float],
    condition_names: list[str],
    title: str,
) -> None:
    """Write the genotype x environment panels.

    Each panel shows one gene's partial-dependence curve within each condition, centred so
    that only the *shape* differences are visible -- a pure condition main effect would
    otherwise dominate the plot and hide the interaction. Curves are the mean over seeds,
    with the seed standard deviation shaded, so the drawn lines match the H^2 in the title.
    """
    # rainbow sampled at the number of conditions, so 18-39 lines stay separable
    colours = plt.cm.rainbow(np.linspace(0, 1, len(condition_names)))
    per_page = 4
    pages = max(1, -(-len(genes) // per_page))

    with PdfPages(path) as pdf:
        for start in range(0, len(genes), per_page):
            chunk = genes[start : start + per_page]
            figure, axes = plt.subplots(2, 2, figsize=(13, 8.5))

            for offset, (axis, gene) in enumerate(zip(axes.ravel(), chunk)):
                plot_gxe_panel(
                    axis, gene, grids[gene], tables[gene], spreads[gene], scores[gene],
                    condition_names, colours, index=start + offset + 1,
                )

            for axis in axes.ravel()[len(chunk) :]:
                axis.axis("off")

            handles, labels = axes.ravel()[0].get_legend_handles_labels()
            figure.legend(
                handles, labels, loc="lower center", ncol=min(8, len(labels)), fontsize=6,
                frameon=False,
            )
            figure.suptitle(
                f"{title}  —  panels {start + 1}-{start + len(chunk)} of {len(genes)} "
                f"(page {start // per_page + 1} of {pages})",
                fontsize=11,
            )
            figure.tight_layout(rect=(0, 0.10, 1, 0.96))
            pdf.savefig(figure)
            plt.close(figure)


def safe_name(name: str) -> str:
    """Filesystem-safe form of a feature name.

    The count encoding names columns ``<GENE>_variant_score=2``, so the ``=`` has to go.
    """
    return re.sub(r"[^0-9A-Za-z._-]", "_", name)


def write_pdp_panels(
    out_dir: Path,
    names: list[str],
    grids: dict[str, ndarray],
    curves: dict[str, ndarray],
    ice: dict[str, ndarray],
) -> None:
    """Write each partial-dependence panel to its own SVG.

    The multi-page PDF is for reading; these are for lifting a single feature into a figure
    without cropping a page. Numbering matches the panel numbers in the PDF.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    for index, name in enumerate(names, start=1):
        figure, axis = plt.subplots(figsize=(5.5, 4.0))
        plot_pdp_page(axis, grids[name], curves[name], ice.get(name), name, index=index)
        figure.tight_layout()
        figure.savefig(out_dir / f"{index:03d}_{safe_name(name)}.svg", format="svg")
        plt.close(figure)


def write_gxe_panels(
    out_dir: Path,
    genes: list[str],
    grids: dict[str, ndarray],
    tables: dict[str, ndarray],
    spreads: dict[str, ndarray],
    scores: dict[str, float],
    condition_names: list[str],
) -> None:
    """Write each genotype x environment panel to its own SVG, with its own legend.

    Unlike the PDF, where one figure-level legend serves four panels, each file has to carry
    the condition key itself to stand alone.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    colours = plt.cm.rainbow(np.linspace(0, 1, len(condition_names)))

    for index, gene in enumerate(genes, start=1):
        figure, axis = plt.subplots(figsize=(7.0, 5.6))
        plot_gxe_panel(
            axis, gene, grids[gene], tables[gene], spreads[gene], scores[gene],
            condition_names, colours, index=index,
        )
        handles, labels = axis.get_legend_handles_labels()
        figure.legend(
            handles, labels, loc="lower center", ncol=min(5, len(labels)), fontsize=5,
            frameon=False,
        )
        figure.tight_layout(rect=(0, 0.16, 1, 1))
        figure.savefig(out_dir / f"{index:03d}_{safe_name(gene)}.svg", format="svg")
        plt.close(figure)


def write_epistasis_heatmap(
    path: Path, genes: list[str], matrix: ndarray, title: str
) -> None:
    """Write the gene x gene H^2 heatmap."""
    size = max(6.0, 0.32 * len(genes) + 3.0)
    figure, axis = plt.subplots(figsize=(size, size * 0.9))

    image = axis.imshow(matrix, cmap="magma", vmin=0.0, origin="upper")
    axis.set_xticks(range(len(genes)), genes, rotation=90, fontsize=6)
    axis.set_yticks(range(len(genes)), genes, fontsize=6)
    axis.set_title(title, fontsize=10)
    figure.colorbar(image, ax=axis, shrink=0.75, label="Friedman's H²")
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


# --------------------------------------------------------------------------------------
# Per-dataset driver
# --------------------------------------------------------------------------------------


def analyse_dataset(
    conf: DictConfig, name: str, model_paths: list[Path], data_path: str, out_dir: Path
) -> None:
    """Run the whole analysis for one dataset across its five seeds.

    Parameters
    ----------
    conf : DictConfig
        The composed configuration.
    name : str
        Dataset name, e.g. ``Bloom2013``.
    model_paths : list[Path]
        Paths to the seeds' ``Boosting.pkl``.
    data_path : str
        Path to the matching feather file.
    out_dir : Path
        Directory to write this dataset's output into.
    """
    settings = conf.analysis
    rng = np.random.default_rng(settings.seed)
    out_dir.mkdir(parents=True, exist_ok=True)

    console.log(f"[{name}] loading {data_path}")
    X, y, conditions, feature_names = load_dataset(data_path)
    console.log(
        f"[{name}] {X.shape[0]} rows x {X.shape[1]} features "
        f"(pdp backend: {settings.pdp_backend})"
    )

    unique_conditions = np.unique(conditions)
    latent_columns = [i for i, n in enumerate(feature_names) if not is_gene(n)]
    unique_latent = np.unique(X[:, latent_columns], axis=0)
    console.log(
        f"[{name}] {len(unique_conditions)} conditions -> "
        f"{len(unique_latent)} distinct chemistry vectors"
    )

    # Background rows, shared by every seed so curves are comparable across seeds.
    sample_idx = stratified_sample(conditions, settings.n_sample, rng)
    sample_conditions = conditions[sample_idx]
    condition_counts = Counter(sample_conditions.tolist())
    kept_conditions = [
        c for c in unique_conditions
        if condition_counts[c] >= settings.min_rows_per_condition
    ]
    if len(kept_conditions) < 2:
        # n_sample // n_conditions fell below the threshold, so the filter would drop
        # everything and silently disable the whole G x E analysis. Keep every condition
        # that is present and say so, rather than emitting an empty result.
        console.log(
            f"[{name}] [yellow]warning[/yellow]: min_rows_per_condition="
            f"{settings.min_rows_per_condition} would keep {len(kept_conditions)} of "
            f"{len(unique_conditions)} conditions at n_sample={settings.n_sample}; "
            "falling back to all conditions. Raise n_sample for denser curves."
        )
        kept_conditions = list(unique_conditions)
    kept_conditions = sorted(kept_conditions)
    condition_codes = np.searchsorted(
        np.asarray(kept_conditions), sample_conditions
    )
    keep_mask = np.isin(sample_conditions, kept_conditions)
    console.log(
        f"[{name}] background sample: {len(sample_idx)} rows across "
        f"{len(kept_conditions)} conditions "
        f"(min {min(condition_counts[c] for c in kept_conditions)} rows/condition)"
    )

    X_sample_raw = X[sample_idx]
    y_sample = y[sample_idx]

    # ---- rank features -------------------------------------------------------------
    gain = np.zeros(X.shape[1])
    models, scalers, aucs = [], [], []

    for model_path in model_paths:
        seed = seed_from_run_dir(Path(model_path))
        model = get_model(Path(model_path))
        scaler, auc = fit_run_scaler(
            X, y, seed, conf.testing.test_frac, settings.verify_auc, model
        )
        gain += model.feature_importance("gain")
        models.append(model)
        scalers.append(scaler)
        if auc is not None:
            aucs.append(auc)
        console.log(
            f"[{name}] seed {seed}: trees={model.num_trees()}"
            + (f", held-out AUC={auc:.4f}" if auc is not None else "")
        )

    if aucs:
        console.log(
            f"[{name}] held-out AUC across seeds: "
            f"{np.mean(aucs):.4f} ± {np.std(aucs):.4f}"
        )

    # Nominate candidates from both cheap rankings. Gain scores how the model was fitted
    # and mean |SHAP| scores it observationally; a feature that either one buries is a
    # feature permutation importance never gets to see. Neither is used to rank --
    # permutation importance still decides the top n_top.
    gain_pool = np.argsort(-gain)[: settings.n_candidates]
    shap_order, shap_values = load_shap_ranking(
        conf.get("shap_path"),
        name,
        conf.model_load_keys.model_type,
        feature_names,
    )

    if shap_order is not None and settings.shap_candidates > 0:
        shap_pool = shap_order[: settings.shap_candidates]
        candidates = np.union1d(gain_pool, shap_pool)
        shap_only = len(np.setdiff1d(shap_pool, gain_pool))
        console.log(
            f"[{name}] candidate pool: {len(gain_pool)} by gain u "
            f"{len(shap_pool)} by |SHAP| = {len(candidates)} "
            f"({len(gain_pool) + len(shap_pool) - len(candidates)} shared, "
            f"{shap_only} nominated by SHAP alone)"
        )
    else:
        candidates = gain_pool
        shap_pool = np.array([], dtype=int)
        console.log(f"[{name}] candidate pool: {len(candidates)} by gain only")

    console.log(f"[{name}] permutation-scoring {len(candidates)} candidates")

    drops = np.zeros(len(candidates))
    for model, scaler in zip(models, scalers):
        X_scaled = scaler.transform(X_sample_raw)
        drops += permutation_importance(
            model, X_scaled, y_sample, candidates, settings.n_permute, rng
        )
        del X_scaled
    drops /= len(models)

    order = np.argsort(-drops)
    ranked = candidates[order]

    # Which cheap ranking nominated each feature. Reported so the effect of the union is
    # auditable rather than buried in the pool size.
    in_gain = set(gain_pool.tolist())
    in_shap = set(shap_pool.tolist())
    source = [
        "both" if c in in_gain and c in in_shap else "gain" if c in in_gain else "shap"
        for c in ranked
    ]

    if shap_values is not None:
        # Rank over every feature, not just the pool, so the number stays comparable
        # across datasets and across changes to n_candidates.
        shap_rank = np.empty(len(shap_values), dtype=int)
        shap_rank[shap_order] = np.arange(1, len(shap_values) + 1)
        shap_column = shap_values[ranked]
        rank_column = shap_rank[ranked]
    else:
        shap_column = [None] * len(ranked)
        rank_column = [None] * len(ranked)

    ranking = pl.DataFrame(
        {
            "Feature": [feature_names[c] for c in ranked],
            "PermutationImportance": drops[order],
            "Gain": gain[ranked],
            "ShapValue": shap_column,
            "ShapRank": rank_column,
            "Source": source,
            "Block": [
                "genotype" if is_gene(feature_names[c]) else "chemistry" for c in ranked
            ],
        }
    )
    ranking.write_csv(out_dir / "feature_ranking.csv")

    ranked_columns = ranked
    top_columns = ranked_columns[: settings.n_top]
    top_names = [feature_names[c] for c in top_columns]
    console.log(
        f"[{name}] top {len(top_names)}: "
        f"{sum(is_gene(n) for n in top_names)} genotype, "
        f"{sum(not is_gene(n) for n in top_names)} chemistry"
    )

    # The genotype block is sparse (~78% zeros) and individually weak, so a single global
    # ranking comes out chemistry-dominated and would leave the G x E analysis -- the whole
    # point of this script -- with a handful of genes. Genes for the interaction sections
    # are therefore ranked within the genotype block, independently of the top-N cut used
    # for the PDP panels.
    gene_columns = [c for c in ranked_columns if is_gene(feature_names[c])][
        : settings.n_gxe_genes
    ]
    gene_names = [feature_names[c] for c in gene_columns]
    console.log(f"[{name}] {len(gene_names)} genes selected for the interaction analysis")

    # Sweep the union once; the PDP pdf plots only the top-N, the interaction sections use
    # the genes, and neither pays for the other.
    sweep_columns = list(top_columns) + [c for c in gene_columns if c not in set(top_columns)]
    sweep_names = [feature_names[c] for c in sweep_columns]

    # ---- 1-D partial dependence, and G x E in the same sweep ------------------------
    # Grid on the background sample, not the full dataset. sklearn derives its grid from
    # the X it is handed, so gridding on full-data support would desynchronise the two and
    # break the equivalence with the builtin backend. The only points this drops are
    # genotype levels absent from the sample -- 5-9 features per dataset, each a level
    # carried by 9-39 of 18k-87k rows, where a partial-dependence estimate has essentially
    # no support anyway.
    grids = {
        n: build_grid(X_sample_raw[:, c], settings.max_grid)
        for n, c in zip(sweep_names, sweep_columns)
    }
    curves: dict[str, ndarray] = {n: np.empty((len(models), len(grids[n]))) for n in sweep_names}
    ice: dict[str, ndarray] = {}
    gxe_tables: dict[str, list[ndarray]] = {}
    gxe_scores: dict[str, list[float]] = {n: [] for n in sweep_names}

    for seed_pos, (model, scaler) in enumerate(zip(models, scalers)):
        X_scaled = scaler.transform(X_sample_raw)

        for name_i, column in zip(sweep_names, sweep_columns):
            grid = grids[name_i]
            grid_scaled = (grid - scaler.mean_[column]) / scaler.scale_[column]
            surface = pdp_sweep(
                model, X_scaled, column, grid_scaled,
                max_grid=settings.max_grid, backend=settings.pdp_backend,
            )

            curves[name_i][seed_pos] = surface.mean(axis=1)

            # Background rows are shared across seeds, so ICE curves can be averaged
            # row-wise rather than taken from a single model.
            chunk = surface[:, : settings.n_ice].T
            ice[name_i] = chunk if seed_pos == 0 else ice[name_i] + chunk

            if name_i in set(gene_names) and len(kept_conditions) > 1:
                table = np.empty((len(grid), len(kept_conditions)))
                for code in range(len(kept_conditions)):
                    member = keep_mask & (condition_codes == code)
                    table[:, code] = surface[:, member].mean(axis=1)
                gxe_tables.setdefault(name_i, []).append(table)
                gxe_scores[name_i].append(
                    gxe_interaction(
                        surface[:, keep_mask],
                        condition_codes[keep_mask],
                        len(kept_conditions),
                    )
                )

        del X_scaled
        console.log(f"[{name}] PDP sweep done for seed position {seed_pos + 1}/{len(models)}")

    for key in ice:
        ice[key] = ice[key] / len(models)

    write_pdp_pdf(
        out_dir / f"pdp_top{len(top_names)}.pdf", top_names, grids, curves, ice,
        f"{name} — partial dependence, top {len(top_names)} features "
        f"(mean ± sd over {len(models)} seeds)",
    )
    write_pdp_panels(out_dir / "pdp_panels", top_names, grids, curves, ice)

    gxe_genes = [n for n in gene_names if n in gxe_tables]
    gxe_mean = {g: float(np.nanmean(gxe_scores[g])) for g in gxe_genes}
    gxe_curves = {g: np.mean(np.stack(gxe_tables[g]), axis=0) for g in gxe_genes}
    gxe_spread = {g: np.std(np.stack(gxe_tables[g]), axis=0) for g in gxe_genes}
    if gxe_genes:
        gxe_genes.sort(key=lambda g: -gxe_mean[g])
        condition_labels = [str(c) for c in kept_conditions]
        write_gxe_pdf(
            out_dir / "gxe_condition_pdp.pdf", gxe_genes, grids, gxe_curves, gxe_spread,
            gxe_mean, condition_labels,
            f"{name} — genotype × environment "
            f"({len(kept_conditions)} conditions, mean ± sd over {len(models)} seeds)",
        )
        write_gxe_panels(
            out_dir / "gxe_panels", gxe_genes, grids, gxe_curves, gxe_spread, gxe_mean,
            condition_labels,
        )

    # ---- genotype x genotype -------------------------------------------------------
    # Every seed contributes, so the epistasis numbers carry the same 5-model support as
    # the marginals and the G x E statistic, and a seed spread can be reported.
    epistasis_names = gene_names[: settings.n_epistasis]
    epistasis_columns = gene_columns[: settings.n_epistasis]
    matrix = np.full((len(epistasis_names), len(epistasis_names)), np.nan)
    spread_matrix = np.full((len(epistasis_names), len(epistasis_names)), np.nan)
    pair_rows = []

    want_ale = settings.interaction_method in ("ale", "both")

    if len(epistasis_names) > 1:
        pair_scores: dict[tuple[str, str], list[float]] = {}
        pair_extra: dict[tuple[str, str], dict[str, int]] = {}

        for seed_pos, (model, scaler) in enumerate(zip(models, scalers)):
            X_scaled = scaler.transform(X_sample_raw)

            if want_ale:
                estimator = BoosterEstimator(model).fit(X_scaled)
                frame = pd.DataFrame(X_scaled, columns=feature_names)
                ale_1d_cache = {n: ale_1d(frame, estimator, n) for n in epistasis_names}

            scaled_grids = {
                n: (grids[n] - scaler.mean_[c]) / scaler.scale_[c]
                for n, c in zip(epistasis_names, epistasis_columns)
            }
            marginals = {
                n: pdp_sweep(
                    model, X_scaled, c, scaled_grids[n],
                    max_grid=settings.max_grid, backend=settings.pdp_backend,
                ).mean(axis=1)
                for n, c in zip(epistasis_names, epistasis_columns)
            }

            for i in range(len(epistasis_names)):
                for j in range(i + 1, len(epistasis_names)):
                    name_a, name_b = epistasis_names[i], epistasis_names[j]
                    key = (name_a, name_b)
                    joint = pdp_2d(
                        model, X_scaled, epistasis_columns[i], epistasis_columns[j],
                        scaled_grids[name_a], scaled_grids[name_b],
                        max_grid=settings.max_grid, backend=settings.pdp_backend,
                    )
                    pair_scores.setdefault(key, []).append(
                        friedman_h2(joint, marginals[name_a], marginals[name_b])
                    )

                    if want_ale and seed_pos == 0:
                        surface = ale_2d(
                            frame, estimator, name_a, name_b, settings.ale_grid_size
                        )
                        col_a, col_b = epistasis_columns[i], epistasis_columns[j]
                        weights = joint_weights(
                            X_sample_raw[:, col_a], X_sample_raw[:, col_b],
                            np.asarray(surface.index, dtype=float) * scaler.scale_[col_a]
                            + scaler.mean_[col_a],
                            np.asarray(surface.columns, dtype=float) * scaler.scale_[col_b]
                            + scaler.mean_[col_b],
                        )
                        pair_extra[key] = {
                            "H2_ALE": h2_ale(
                                surface, ale_1d_cache[name_a], ale_1d_cache[name_b], weights
                            ),
                            "EmptyCells": int((weights == 0).sum()),
                            "Cells": int(weights.size),
                        }

            del X_scaled
            console.log(
                f"[{name}] epistasis done for seed position "
                f"{seed_pos + 1}/{len(models)}"
            )

        for i in range(len(epistasis_names)):
            for j in range(i + 1, len(epistasis_names)):
                key = (epistasis_names[i], epistasis_names[j])
                values = np.asarray(pair_scores[key], dtype=float)
                mean = float(np.nanmean(values))
                sd = float(np.nanstd(values))
                matrix[i, j] = matrix[j, i] = mean
                spread_matrix[i, j] = spread_matrix[j, i] = sd
                record = {"FeatureA": key[0], "FeatureB": key[1], "H2": mean, "H2_sd": sd}
                record.update(pair_extra.get(key, {}))
                pair_rows.append(record)

        write_epistasis_heatmap(
            out_dir / "epistasis_heatmap.png", epistasis_names, np.nan_to_num(matrix),
            f"{name} — gene × gene Friedman's H² "
            f"(top {len(epistasis_names)} genes, mean of {len(models)} seeds)",
        )

    # ---- one CSV holding every interaction statistic --------------------------------
    records = [
        {
            "Kind": "gene_x_condition",
            "FeatureA": gene,
            "FeatureB": "Condition",
            "H2": gxe_mean[gene],
            "H2_sd": float(np.nanstd(gxe_scores[gene])),
        }
        for gene in gxe_genes
    ]
    records += [dict(r, Kind="gene_x_gene") for r in pair_rows]

    if records:
        (
            pl.DataFrame(records)
            .sort("H2", descending=True, nulls_last=True)
            .write_csv(out_dir / "interaction_strength.csv")
        )

    # Conditions that share a chemistry vector are indistinguishable to the model.
    collisions = len(unique_conditions) - len(unique_latent)
    pl.DataFrame(
        {
            "Dataset": [name],
            "Rows": [X.shape[0]],
            "Features": [X.shape[1]],
            "Conditions": [len(unique_conditions)],
            "DistinctChemistryVectors": [len(unique_latent)],
            "CollidingConditions": [collisions],
            "ConditionsUsed": [len(kept_conditions)],
            "HeldOutAUCMean": [float(np.mean(aucs)) if aucs else None],
            "HeldOutAUCSd": [float(np.std(aucs)) if aucs else None],
            "MedianGxE": [float(np.nanmedian(list(gxe_mean.values()))) if gxe_mean else None],
            "MedianGxG": [
                float(np.nanmedian([r["H2"] for r in pair_rows])) if pair_rows else None
            ],
        }
    ).write_csv(out_dir / "summary.csv")

    if collisions:
        console.log(
            f"[{name}] [yellow]warning[/yellow]: {collisions} condition(s) share a "
            "chemistry vector with another condition"
        )

    del X, y
    console.log(f"[{name}] [green]done[/green] -> {out_dir}")


@hydra.main(  # pyrefly: ignore
    config_path="../configs/", version_base="1.3", config_name="interaction"
)
def main(conf: DictConfig) -> None:
    """Entry point. Runs :func:`analyse_dataset` for every configured dataset."""
    out_path = Path(conf.out_path)
    out_path.mkdir(parents=True, exist_ok=True)

    model_paths = get_model_paths(**conf.model_load_keys)

    for name, paths in model_paths.items():
        if not paths:
            console.log(f"[{name}] [red]no models found[/red], skipping")
            continue
        analyse_dataset(
            conf, name, sorted(paths), conf.data_paths.get(name), out_path / name
        )


if __name__ == "__main__":
    main()
