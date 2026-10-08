"""

``src/shap_analysis.py`` explains a model over its whole dataset at once, so what comes out
is a gene ranking averaged across every drug and carbon source the panel was grown in, with
the 256 ``latent_*`` chemical-embedding columns taking 68-80% of the total attribution. This
script explains one condition at a time instead.
"""

from pathlib import Path

import hydra
import numpy as np
import polars as pl
import shap
from dotenv import load_dotenv
from numpy import ndarray
from omegaconf import DictConfig
from rich.console import Console
from sklearn.preprocessing import StandardScaler
from utils import get_model, get_model_paths

load_dotenv()

console = Console()

META_COLUMNS = ["Condition", "Strain", "Phenotype"]
LATENT_PREFIX = "latent"

# Default ceiling on latent attribution, overridable as `latent_tolerance`. Exact
# interventional TreeSHAP returns 0.000e+00 for the latent block on every seed measured, so
# anything above floating-point noise means the premise of the analysis has failed and the
# run should stop rather than write a misleading table. Saabas does NOT hold to this -- it
# leaks 4e-2 to 9e-2 on four of five Bloom2013 seeds, and a larger background does not fix
# it -- so the companion run raises the ceiling deliberately and logs what leaked.
DEFAULT_LATENT_TOLERANCE = 1e-10


def load_matrix(data_path: str) -> tuple[ndarray, ndarray, list[str], list[str]]:
    """Load one dataset and standardise it the way the models were trained.

    The scaler is fitted across every row, matching ``shap_analysis.py``. Subsetting happens
    afterwards, on the already-standardised matrix.

    Parameters
    ----------
    data_path : str
        Path to the feather the models were trained on.

    Returns
    -------
    tuple[ndarray, ndarray, list[str], list[str]]
        The standardised feature matrix, the ``Condition`` column, all feature names in
        column order, and the gene feature names alone.
    """
    data = pl.read_ipc(data_path)
    features = data.drop(META_COLUMNS)
    names = features.columns

    matrix = StandardScaler().fit_transform(features.to_numpy().astype(np.float64))
    genes = [name for name in names if not name.startswith(LATENT_PREFIX)]

    return matrix, data["Condition"].to_numpy(), names, genes


def fold_shap(
    conf: DictConfig,
    model,
    matrix: ndarray,
    conditions: ndarray,
    order: list[str],
    gene_index: ndarray,
    latent_index: ndarray,
) -> ndarray:
    """Mean absolute SHAP per gene, per condition, for a single fold's model.

    Parameters
    ----------
    conf : DictConfig
        Composed configuration; supplies ``explainer``, ``n_background``,
        ``background_seed``, ``approximate`` and ``check_additivity``.
    model : lgb.Booster
        The fold's fitted model.
    matrix : ndarray
        Standardised feature matrix for the whole dataset.
    conditions : ndarray
        The dataset's ``Condition`` column, aligned to ``matrix`` rows.
    order : list[str]
        Conditions to explain, in output order.
    gene_index : ndarray
        Column indices of the gene features.
    latent_index : ndarray
        Column indices of the ``latent_*`` features, checked against
        ``conf.latent_tolerance``.

    Returns
    -------
    ndarray
        Array of shape ``(len(order), len(gene_index))``.

    Raises
    ------
    ValueError
        If a condition's latent SHAP is not zero, which means its chemical embedding varies
        within the condition and the explanation cannot be read as gene-only.
    """
    rng = np.random.default_rng(conf.background_seed)
    out = np.zeros((len(order), gene_index.size))

    for position, condition in enumerate(order):
        rows = matrix[conditions == condition]
        size = min(int(conf.n_background), rows.shape[0])
        background = rows[rng.choice(rows.shape[0], size=size, replace=False)]

        explainer: shap.TreeExplainer = hydra.utils.instantiate(
            conf.explainer, model, background
        )
        values = explainer.shap_values(
            rows, approximate=conf.approximate, check_additivity=conf.check_additivity
        )

        tolerance = float(conf.get("latent_tolerance", DEFAULT_LATENT_TOLERANCE))
        leaked = np.abs(values[:, latent_index]).max()
        if leaked > tolerance:
            raise ValueError(
                f"condition {condition!r}: latent features carry SHAP up to {leaked:.3e}, "
                f"above the {tolerance:.1e} ceiling. Either the chemical embedding is not "
                "constant within this condition, or the estimator does not zero it -- "
                "Saabas does not, and needs latent_tolerance raised explicitly"
            )

        out[position] = np.abs(values[:, gene_index]).mean(axis=0)
        console.log(
            f"  {condition:<24} {rows.shape[0]:>5} rows, background {size:>4}, "
            f"max |SHAP| {out[position].max():.4f}, latent leak {leaked:.2e}"
        )

    return out


def analyse_dataset(conf: DictConfig, model_paths: list[Path], data_path: str) -> pl.DataFrame:
    """Explain every fold of one dataset and average the folds.

    Parameters
    ----------
    conf : DictConfig
        Composed configuration.
    model_paths : list[Path]
        The dataset's ``Boosting.pkl`` paths, one per seed.
    data_path : str
        Path to the dataset's feather.

    Returns
    -------
    pl.DataFrame
        Long frame of ``Condition``, ``Feature``, ``Value``, ``Rows``.
    """
    matrix, conditions, names, genes = load_matrix(data_path)

    order = list(conf.conditions) if conf.get("conditions") else sorted(set(conditions))
    missing = [c for c in order if c not in set(conditions)]
    if missing:
        raise ValueError(f"conditions absent from {data_path}: {missing}")

    gene_index = np.array([i for i, name in enumerate(names) if not name.startswith(LATENT_PREFIX)])
    latent_index = np.array([i for i, name in enumerate(names) if name.startswith(LATENT_PREFIX)])
    console.log(
        f"{matrix.shape[0]} rows, {gene_index.size} genes, {latent_index.size} latent, "
        f"{len(order)} conditions, {len(model_paths)} folds"
    )

    total = np.zeros((len(order), gene_index.size))
    for fold, model_path in enumerate(model_paths):
        console.log(f"fold {fold}: {Path(model_path).parent.name}")
        total += fold_shap(
            conf,
            get_model(model_path),
            matrix,
            conditions,
            order,
            gene_index,
            latent_index,
        )
    total /= len(model_paths)

    rows = {condition: int((conditions == condition).sum()) for condition in order}

    return pl.DataFrame(
        {
            "Condition": np.repeat(order, gene_index.size),
            "Feature": np.tile(genes, len(order)),
            "Value": total.ravel(),
            "Rows": np.repeat([rows[c] for c in order], gene_index.size),
        }
    ).with_columns(pl.col("Condition").cast(pl.Categorical))


@hydra.main(  # pyrefly: ignore
    config_path="../configs/", version_base="1.3", config_name="interpret_conditionwise"
)
def main(conf: DictConfig) -> None:
    """Entry point."""
    out_path = Path(conf.out_path)
    out_path.mkdir(parents=True, exist_ok=True)
    model_type = conf.model_load_keys.model_type
    model_paths = get_model_paths(**conf.model_load_keys)

    for model_name, paths in model_paths.items():
        paths = sorted(paths)
        if not paths:
            console.log(f"[{model_name}] [red]no models found[/red], skipping")
            continue

        console.log(f"[bold]{model_name}[/bold]")
        frame = analyse_dataset(conf, paths, conf.data_paths.get(model_name))

        dest = out_path / f"{model_name}_{model_type}_conditionwise_shap.parquet"
        frame.write_parquet(dest)
        console.log(f"[{model_name}] -> {dest} ({len(frame)} rows)")

    console.log(f"[green]done[/green] -> {out_path}")


if __name__ == "__main__":
    main()
