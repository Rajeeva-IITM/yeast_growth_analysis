"""Distribution of gene-level encodings across the gene x strain matrix.
"""

import argparse
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
from dotenv import load_dotenv  # noqa: E402

load_dotenv()

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.size"] = 10

LEVELS = (0, 1, 2, 3)

# The max encoding lives in the unsuffixed file; _sum and _count are the other encodings.
SOURCES = {
    "Bloom2013": "bloom2013.parquet",
    "Bloom2015": "bloom2015_classification.parquet",
    "Bloom2019_BYxRM": "bloom2019_BYxRM.parquet",
}

COLOURS = {"Bloom2013": "#0072B2", "Bloom2015": "#009E73", "Bloom2019_BYxRM": "#D55E00"}


def load_matrix(path: Path) -> np.ndarray:
    """One row per strain, one column per gene.

    Parameters
    ----------
    path : Path
        A per-strain dataset parquet.

    Returns
    -------
    np.ndarray
        The ``(n_strains, n_genes)`` encoding matrix.
    """
    genes = [
        name
        for name in pq.read_schema(path).names
        if name.startswith("Y") and "latent" not in name
    ]
    frame = (
        pl.scan_parquet(path)
        .unique(subset=["Strain"], keep="first")
        .select(genes)
        .collect(streaming=True)
    )
    return frame.to_numpy()


def summarise(matrix: np.ndarray, dataset: str) -> list[dict]:
    """Counts of each encoding level, over all genes and over segregating genes.

    Parameters
    ----------
    matrix : np.ndarray
        Encoding matrix from :func:`load_matrix`.
    dataset : str
        Dataset name, copied into every row.

    Returns
    -------
    list[dict]
        One row per (scope, level), plus per-gene and per-strain summaries.
    """
    segregating = matrix.max(axis=0) > 0
    rows = []

    for scope, block in (
        ("all genes", matrix),
        ("segregating genes", matrix[:, segregating]),
    ):
        counts = np.bincount(block.ravel(), minlength=len(LEVELS))[: len(LEVELS)]
        total = counts.sum()
        for level in LEVELS:
            rows.append(
                {
                    "Dataset": dataset,
                    "Scope": scope,
                    "Encoding": level,
                    "Cells": int(counts[level]),
                    "Percent": 100.0 * counts[level] / total,
                    "GenesAtMax": int((matrix.max(axis=0) == level).sum()),
                    "MeanGenesPerStrain": float(
                        np.mean((matrix == level).sum(axis=1))
                    ),
                }
            )

    return rows


def draw(frame: pl.DataFrame, out_path: Path) -> None:
    """Grouped bars, one panel per scope, log y because level 0 dwarfs level 3.

    Parameters
    ----------
    frame : pl.DataFrame
        Output of :func:`summarise` for every dataset.
    out_path : Path
        Destination ``.svg``.

    Returns
    -------
    None
    """
    scopes = ["all genes", "segregating genes"]
    figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.4), sharey=True)
    width = 0.26

    for axis, scope in zip(axes, scopes):
        block = frame.filter(pl.col("Scope") == scope)
        for offset, dataset in enumerate(SOURCES):
            values = (
                block.filter(pl.col("Dataset") == dataset).sort("Encoding")["Percent"].to_list()
            )
            positions = np.arange(len(LEVELS)) + (offset - 1) * width
            bars = axis.bar(
                positions, values, width, color=COLOURS[dataset], label=dataset
            )
            for bar, value in zip(bars, values):
                axis.text(
                    bar.get_x() + bar.get_width() / 2,
                    value * 1.12,
                    f"{value:.2f}" if value < 10 else f"{value:.1f}",
                    ha="center",
                    fontsize=6.5,
                )

        axis.set_yscale("log")
        axis.set_xticks(range(len(LEVELS)), [str(level) for level in LEVELS])
        axis.set_xlabel("gene-level encoding")
        axis.set_title(scope, fontsize=10)
        axis.grid(axis="y", alpha=0.2, linewidth=0.5)
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)

    axes[0].set_ylabel("% of gene × strain entries")
    axes[0].legend(frameon=False, fontsize=8, loc="lower left")
    figure.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    """Entry point."""
    project_dir = Path(os.environ["PROJECT_DIR"])

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path.home() / "projects/bloom_data/datasets",
        help="directory holding the per-strain dataset parquets",
    )
    parser.add_argument(
        "--out-dir", type=Path, default=project_dir / "results/encoding_distribution"
    )
    args = parser.parse_args()

    rows = []
    for dataset, filename in SOURCES.items():
        matrix = load_matrix(args.data_dir / filename)
        rows.extend(summarise(matrix, dataset))
        print(
            f"{dataset}: {matrix.shape[0]} strains x {matrix.shape[1]} genes, "
            f"{int((matrix.max(axis=0) > 0).sum())} segregating"
        )
        del matrix

    frame = pl.DataFrame(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    frame.write_csv(args.out_dir / "encoding_distribution.csv")
    draw(frame, args.out_dir / "encoding_distribution.svg")
    print(f"-> {args.out_dir}")


if __name__ == "__main__":
    main()
