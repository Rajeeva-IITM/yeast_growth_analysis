"""Held-out predictions and per-condition metrics for the trained models.

"""

from pathlib import Path

import hydra
import numpy as np
import polars as pl
from dotenv import load_dotenv
from numpy import ndarray
from omegaconf import DictConfig
from rich.console import Console
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from utils import get_data, get_model, get_model_paths

load_dotenv()

console = Console()


def seed_from_run_dir(model_path: Path) -> int:
    """Recover the training seed from a run directory named ``<Prefix>_<Name>_<seed>_<suffix>``.

    Mirrors the function of the same name in ``src/interaction_analysis.py``; the two must
    agree, since both reproduce the same training split.

    Parameters
    ----------
    model_path : Path
        Path to ``Boosting.pkl`` inside the run directory.

    Returns
    -------
    int
        The seed encoded in the directory name.
    """
    for token in reversed(model_path.parent.name.split("_")):
        if token.isdigit():
            return int(token)
    raise ValueError(f"no seed found in run directory name {model_path.parent.name!r}")


def run_split(n_rows: int, seed: int, test_frac: float) -> tuple[ndarray, ndarray]:
    """Reproduce one run's train/test partition.

    ``tune_model.py`` calls ``train_test_split(test_size=test_frac, random_state=seed)``.
    Splitting an index array of the same length reproduces the identical partition without
    copying a multi-GB feature matrix.

    Parameters
    ----------
    n_rows : int
        Number of rows in the training dataset.
    seed : int
        The run's seed.
    test_frac : float
        Test fraction, from ``conf.testing.test_frac``. Must match ``conf.yaml``.

    Returns
    -------
    tuple[ndarray, ndarray]
        Train and test row indices.
    """
    return train_test_split(
        np.arange(n_rows), test_size=test_frac, random_state=seed
    )


def source_fits(
    data_path: str, run_type: str, seeds: list[int], test_frac: float
) -> dict[int, tuple[StandardScaler, ndarray]]:
    """Fit each seed's scaler on its own training rows, and record its held-out rows.

    Done once per model family so the source matrix is loaded a single time rather than
    once per (model, dataset) pair.

    Parameters
    ----------
    data_path : str
        Path to the feather the models were trained on.
    run_type : str
        ``conf.run_type``; ``full`` for these runs.
    seeds : list[int]
        Training seeds, one per run directory.
    test_frac : float
        Test fraction.

    Returns
    -------
    dict[int, tuple[StandardScaler, ndarray]]
        Per seed, the scaler fitted on that run's training rows and its test indices.
    """
    X, _ = get_data(data_path=data_path, run_type=run_type, return_as_Xy=True)

    fits = {}
    for seed in seeds:
        train_idx, test_idx = run_split(X.shape[0], seed, test_frac)
        fits[seed] = (StandardScaler().fit(X[train_idx]), test_idx)

    del X
    return fits


def predict_frame(
    model_paths: list[Path],
    fits: dict[int, tuple[StandardScaler, ndarray]],
    data_path: str,
    run_type: str,
    matched: bool,
) -> pl.DataFrame:
    """Score every seed of one model family against one dataset.

    Parameters
    ----------
    model_paths : list[Path]
        The family's ``Boosting.pkl`` paths, one per seed.
    fits : dict[int, tuple[StandardScaler, ndarray]]
        Output of :func:`source_fits` for that family.
    data_path : str
        Path to the dataset being scored.
    run_type : str
        ``conf.run_type``.
    matched : bool
        True when the model was trained on this dataset, in which case only that seed's
        held-out rows are scored. False for cross-dataset pairs, where no row was seen in
        training and all of them are used.

    Returns
    -------
    pl.DataFrame
        ``Phenotype``, ``Preds`` (0/1), ``PredProba``, ``Condition``, ``Fold``, ``Split``.
    """
    X, y = get_data(data_path=data_path, run_type=run_type, return_as_Xy=True)
    conditions = pl.read_ipc(data_path, columns=["Condition"])["Condition"].to_numpy()

    frames = []
    for fold, model_path in enumerate(model_paths):
        seed = seed_from_run_dir(Path(model_path))
        scaler, test_idx = fits[seed]
        rows = test_idx if matched else np.arange(X.shape[0])

        proba = get_model(model_path).predict(scaler.transform(X[rows]))
        frames.append(
            pl.DataFrame(
                {
                    "Phenotype": y[rows],
                    "Preds": np.where(proba > 0.5, 1, 0),
                    "PredProba": proba,
                    "Condition": conditions[rows],
                    "Fold": np.full(len(rows), fold, dtype=np.int32),
                    "Split": np.full(
                        len(rows), "holdout" if matched else "cross", dtype=object
                    ),
                }
            )
        )

    del X
    return pl.concat(frames)


def score_group(conf: DictConfig, group: pl.DataFrame) -> dict[str, float]:
    """Apply every configured metric to one Condition x Fold group.

    Metrics named in ``conf.proba_metrics`` are handed the probability; the rest get the
    0/1 label. A group that holds a single class cannot support a threshold-free metric, so
    that metric is recorded as null rather than raising.

    Parameters
    ----------
    conf : DictConfig
        Composed configuration, supplying ``metrics`` and ``proba_metrics``.
    group : pl.DataFrame
        Rows for one condition and fold.

    Returns
    -------
    dict[str, float]
        Metric name to value.
    """
    proba_metrics = set(conf.get("proba_metrics", []))

    scores = {}
    for metric in conf.metrics:
        predicted = group["PredProba"] if metric in proba_metrics else group["Preds"]
        try:
            scores[metric] = float(
                hydra.utils.call(conf.metrics.get(metric), group["Phenotype"], predicted)
            )
        except ValueError:
            scores[metric] = None

    return scores


def eval_model(conf: DictConfig, pred_df: pl.DataFrame) -> pl.DataFrame:
    """Per-condition, per-fold metrics.

    Parameters
    ----------
    conf : DictConfig
        Composed configuration.
    pred_df : pl.DataFrame
        Output of :func:`predict_frame`.

    Returns
    -------
    pl.DataFrame
        One row per condition and fold.
    """
    rows = []
    for (condition, fold), group in pred_df.group_by(
        ["Condition", "Fold"], maintain_order=True
    ):
        rows.append(
            {"Condition": condition, "Fold": fold, "Rows": len(group)}
            | score_group(conf, group)
        )

    return pl.DataFrame(rows)


def get_results(conf: DictConfig) -> None:
    """Score every model family against every dataset and write the output.

    Parameters
    ----------
    conf : DictConfig
        Composed configuration.

    Returns
    -------
    None
    """
    out_path = Path(conf.out_path)
    out_path.mkdir(parents=True, exist_ok=True)
    model_type = conf.model_load_keys.model_type
    model_paths = get_model_paths(**conf.model_load_keys)

    for model_name, paths in model_paths.items():
        paths = sorted(paths)
        if not paths:
            console.log(f"[{model_name}] [red]no models found[/red], skipping")
            continue

        seeds = [seed_from_run_dir(Path(p)) for p in paths]
        console.log(f"[{model_name}] {len(paths)} seeds {seeds}; fitting run scalers")
        fits = source_fits(
            conf.data_paths.get(model_name),
            conf.run_type,
            seeds,
            conf.testing.test_frac,
        )

        for data_name, data_path in conf.data_paths.items():
            matched = data_name == model_name
            pred_df = predict_frame(
                paths, fits, data_path, conf.run_type, matched
            )
            metric_df = eval_model(conf, pred_df)

            pred_df.write_parquet(
                out_path / f"predictions_{model_name}_{data_name}_{model_type}.parquet"
            )
            metric_df.write_csv(
                out_path / f"metrics_{model_name}_{data_name}_{model_type}.csv"
            )

            summary = ", ".join(
                f"{m}={metric_df[m].mean():.4f}"
                for m in conf.metrics
                if metric_df[m].dtype.is_numeric()
            )
            console.log(
                f"[{model_name}] -> {data_name} "
                f"({'held-out' if matched else 'cross'}, {len(pred_df)} rows): {summary}"
            )

        del fits


@hydra.main(  # pyrefly: ignore
    config_path="../configs/", version_base="1.3", config_name="eval_encodings_holdout"
)
def main(conf: DictConfig) -> None:
    """Entry point."""
    get_results(conf)
    console.log(f"[green]done[/green] -> {conf.out_path}")


if __name__ == "__main__":
    main()
