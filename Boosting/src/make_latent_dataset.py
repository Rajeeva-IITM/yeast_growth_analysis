"""Extract the per-condition chemistry embedding into a standalone table.
"""

import argparse
import os
from pathlib import Path

import numpy as np
import polars as pl
from dotenv import load_dotenv

load_dotenv()

N_LATENT = 256
LATENT_COLUMNS = [f"latent_{i}" for i in range(N_LATENT)]

# dataset name -> (feather stem, label column in condition_chemistry.csv)
DATASETS = {
    "Bloom2013": ("bloom2013", "LabelBloom2013"),
    "Bloom2015": ("bloom2015", "LabelBloom2015"),
    "Bloom2019_BYxRM": ("bloom2019_BYxRM", "LabelBloom2019"),
}

# Functional classes, keyed by the canonical `Key` column of condition_chemistry.csv.
# Assignment is by role in the assay, not by formula: SDS and paraquat dichloride are salts
# chemically but act as a detergent and a redox cycler; cisplatin is a metal complex acting
# as a DNA crosslinker; EGTA is a chelator. Sorbitol is Bloom's osmotic stressor rather than
# a carbon source. Every one of these is listed under `other`.
CLASSES = {
    "salt": [
        "CaCl2", "CdCl2", "CoCl2", "CuSO4", "LiCl", "MgCl2", "MgSO4", "MnSO4",
    ],
    "carbon_source": [
        "ethanol", "fructose", "galactose", "glycerol", "lactate", "lactose",
        "maltose", "mannose", "raffinose", "sucrose", "trehalose", "xylose",
    ],
    "other": [
        "4NQO", "6-azauracil", "berbamine", "caffeine", "cisplatin", "congo_red",
        "cycloheximide", "diamide", "egta", "fluconazole", "fluorocytosine",
        "fluorouracil", "formamide", "H2O2", "hydroquinone", "hydroxybenzaldehyde",
        "hydroxyurea", "indoleacetic_acid", "menadione", "methotrexate", "neomycin",
        "paraquat", "SDS", "sorbitol", "tunicamycin", "zeocin",
    ],
}

# Counts asserted below. They come from the data, not from the plan: changing them without
# a corresponding change in the feathers means the merge broke.
EXPECTED_CONDITIONS = {"Bloom2013": 39, "Bloom2015": 18, "Bloom2019_BYxRM": 31}
EXPECTED_CLASS_SIZES = {"salt": 8, "carbon_source": 12, "other": 26}
EXPECTED_KEYS = 46
EXPECTED_DISTINCT_VECTORS = 44


def vector_id(vector: np.ndarray) -> bytes:
    """Hashable identity for one latent vector.

    Rounded at 1e-12 before hashing. The duplicate pairs in these files are bit-identical,
    so the rounding is a guard against a future rebuild introducing float noise, not a
    tolerance the current data needs.

    Parameters
    ----------
    vector : np.ndarray
        One 256-element latent vector.

    Returns
    -------
    bytes
        Byte identity of the rounded vector.
    """
    return np.round(np.asarray(vector, dtype=np.float64), 12).tobytes()


def class_of(key: str) -> str:
    """Functional class for one canonical compound key.

    Parameters
    ----------
    key : str
        Canonical key from ``condition_chemistry.csv``.

    Returns
    -------
    str
        One of ``salt``, ``carbon_source``, ``other``.
    """
    for name, members in CLASSES.items():
        if key in members:
            return name
    raise KeyError(f"no class assigned to compound key {key!r}")


def label_to_key(chemistry: pl.DataFrame) -> dict[str, str]:
    """Map every raw ``Condition`` string used in the feathers to its canonical key.

    Parameters
    ----------
    chemistry : pl.DataFrame
        ``results/condition_chemistry.csv``.

    Returns
    -------
    dict[str, str]
        Raw label to canonical key. Built from the three ``Label*`` columns, which record
        each dataset's own spelling (``copper`` vs ``CuSO4``, ``Diamide`` vs ``diamide``).
    """
    mapping = {}
    for row in chemistry.iter_rows(named=True):
        for _, label_column in DATASETS.values():
            label = row[label_column]
            if label:
                mapping[label] = row["Key"]
    return mapping


def read_conditions(data_dir: Path, stem: str) -> tuple[list[str], np.ndarray]:
    """Read one dataset's distinct conditions and their latent vectors.

    Only the ``Condition`` column and the 256 latent columns are read -- 257 of 6273.

    Parameters
    ----------
    data_dir : Path
        Directory holding the bloom_256 feathers.
    stem : str
        Feather stem, e.g. ``bloom2013``.

    Returns
    -------
    tuple[list[str], np.ndarray]
        Raw condition labels and the matching ``(n_conditions, 256)`` matrix.
    """
    frame = pl.read_ipc(
        data_dir / f"{stem}_max.feather", columns=["Condition"] + LATENT_COLUMNS
    ).unique(subset=["Condition"], maintain_order=True)

    matrix = frame.select(LATENT_COLUMNS).to_numpy().astype(np.float64)
    return frame["Condition"].to_list(), matrix


def collect(data_dir: Path, mapping: dict[str, str]) -> tuple[dict, dict, list[dict]]:
    """Gather every (dataset, label, vector) triple and fold it onto canonical keys.

    Parameters
    ----------
    data_dir : Path
        Directory holding the bloom_256 feathers.
    mapping : dict[str, str]
        Raw label to canonical key, from :func:`label_to_key`.

    Returns
    -------
    tuple[dict, dict, list[dict]]
        Per-key record (vector, labels, membership), per-dataset distinct-vector counts,
        and any rows for the collision audit arising from a key with two vectors.
    """
    records: dict[str, dict] = {}
    per_dataset_vectors: dict[str, int] = {}
    conflicts: list[dict] = []

    for dataset, (stem, _) in DATASETS.items():
        labels, matrix = read_conditions(data_dir, stem)

        unmapped = sorted(set(labels) - set(mapping))
        if unmapped:
            raise KeyError(
                f"{dataset}: condition labels absent from condition_chemistry.csv: {unmapped}"
            )
        if len(labels) != EXPECTED_CONDITIONS[dataset]:
            raise ValueError(
                f"{dataset}: expected {EXPECTED_CONDITIONS[dataset]} conditions, "
                f"read {len(labels)}"
            )

        per_dataset_vectors[dataset] = len({vector_id(row) for row in matrix})

        for label, vector in zip(labels, matrix):
            key = mapping[label]
            record = records.setdefault(
                key, {"vector": None, "source_label": None, "labels": set(), "datasets": set()}
            )
            record["labels"].add(label)
            record["datasets"].add(dataset)

            if record["vector"] is None:
                record["vector"], record["source_label"] = vector, label
            elif vector_id(vector) != vector_id(record["vector"]):
                # One key, two vectors. Prefer the label that matches the key itself.
                kept, dropped = record["vector"], vector
                kept_label, dropped_label = record["source_label"], label
                if label == key:
                    kept, dropped = vector, record["vector"]
                    kept_label, dropped_label = label, record["source_label"]
                    record["vector"], record["source_label"] = vector, label

                distance = float(np.linalg.norm(kept - dropped))
                cosine = float(
                    kept @ dropped / (np.linalg.norm(kept) * np.linalg.norm(dropped))
                )
                conflicts.append(
                    {
                        "Type": "conflicting_vector",
                        "Keys": key,
                        "RawLabels": f"{kept_label}|{dropped_label}",
                        "Detail": (
                            f"kept '{kept_label}', dropped '{dropped_label}'; "
                            f"L2={distance:.6f} cos={cosine:.6f}"
                        ),
                    }
                )

    return records, per_dataset_vectors, conflicts


def build_frame(records: dict, chemistry: pl.DataFrame) -> pl.DataFrame:
    """Assemble the 46-row output table.

    Parameters
    ----------
    records : dict
        Per-key records from :func:`collect`.
    chemistry : pl.DataFrame
        ``results/condition_chemistry.csv``, for ``Name`` / ``CID`` / ``MolecularFormula``.

    Returns
    -------
    pl.DataFrame
        One row per canonical compound, latent columns last.
    """
    meta = {row["Key"]: row for row in chemistry.iter_rows(named=True)}

    ids: dict[bytes, int] = {}
    rows = []
    for key in sorted(records, key=str.lower):
        record = records[key]
        identity = vector_id(record["vector"])
        ids.setdefault(identity, len(ids))

        rows.append(
            {
                "Key": key,
                "Name": meta[key]["Name"],
                "Class": class_of(key),
                "CID": meta[key]["CID"],
                "MolecularFormula": meta[key]["MolecularFormula"],
                "SourceLabel": record["source_label"],
                "RawLabels": ",".join(sorted(record["labels"])),
                "VectorID": ids[identity],
            }
            | {dataset: int(dataset in record["datasets"]) for dataset in DATASETS}
            | {name: float(value) for name, value in zip(LATENT_COLUMNS, record["vector"])}
        )

    order = (
        ["Key", "Name", "Class", "CID", "MolecularFormula", "SourceLabel", "RawLabels", "VectorID"]
        + list(DATASETS)
        + LATENT_COLUMNS
    )
    return pl.DataFrame(rows).select(order)


def audit(frame: pl.DataFrame, conflicts: list[dict]) -> pl.DataFrame:
    """Rows describing every place two compounds or two labels share or contest a vector.

    Parameters
    ----------
    frame : pl.DataFrame
        Output of :func:`build_frame`.
    conflicts : list[dict]
        Conflict rows from :func:`collect`.

    Returns
    -------
    pl.DataFrame
        The collision audit.
    """
    rows = []
    for (identity,), group in frame.group_by(["VectorID"], maintain_order=True):
        if group.height > 1:
            rows.append(
                {
                    "Type": "identical_vector",
                    "Keys": ",".join(sorted(group["Key"].to_list())),
                    "RawLabels": ",".join(sorted(group["RawLabels"].to_list())),
                    "Detail": f"VectorID {identity}: {group.height} compounds share one vector",
                }
            )

    return pl.DataFrame(rows + conflicts)


def cross_check(project_dir: Path, per_dataset_vectors: dict[str, int]) -> None:
    """Compare the distinct-vector counts against those the interaction pipeline recorded.

    ``src/interaction_analysis.py`` writes ``DistinctChemistryVectors`` per dataset. It reads
    the same feathers by a different route, so agreement is an independent check that this
    script's merge is right. Missing summaries are skipped; a mismatch is fatal.

    Parameters
    ----------
    project_dir : Path
        Repository root.
    per_dataset_vectors : dict[str, int]
        Counts derived here.

    Returns
    -------
    None
    """
    base = project_dir / "results/interactions/no_linear_tree/sigma_0.5/max/pdp"
    for dataset, count in per_dataset_vectors.items():
        summary = base / dataset / "summary.csv"
        if not summary.exists():
            print(f"  {dataset}: {count} distinct vectors (no summary.csv to compare)")
            continue

        recorded = int(pl.read_csv(summary)["DistinctChemistryVectors"][0])
        if recorded != count:
            raise ValueError(
                f"{dataset}: derived {count} distinct chemistry vectors but "
                f"{summary} records {recorded}"
            )
        print(f"  {dataset}: {count} distinct vectors, matches interaction_analysis")


def check(frame: pl.DataFrame) -> None:
    """Assert every measured property of the table.

    Parameters
    ----------
    frame : pl.DataFrame
        Output of :func:`build_frame`.

    Returns
    -------
    None
    """
    if frame.height != EXPECTED_KEYS:
        raise ValueError(f"expected {EXPECTED_KEYS} compounds, built {frame.height}")

    latent = [name for name in frame.columns if name.startswith("latent_")]
    if len(latent) != N_LATENT:
        raise ValueError(f"expected {N_LATENT} latent columns, found {len(latent)}")

    for dataset, expected in EXPECTED_CONDITIONS.items():
        # Membership counts canonical compounds, so it is below the raw condition count
        # wherever a dataset spells one compound two ways -- it never is, but assert the
        # relationship rather than assuming it.
        found = int(frame[dataset].sum())
        if found > expected:
            raise ValueError(f"{dataset}: {found} compounds exceeds {expected} conditions")

    sizes = dict(
        frame.group_by("Class").len().iter_rows()  # pyrefly: ignore
    )
    if sizes != EXPECTED_CLASS_SIZES:
        raise ValueError(f"class sizes {sizes} != {EXPECTED_CLASS_SIZES}")

    distinct = frame["VectorID"].n_unique()
    if distinct != EXPECTED_DISTINCT_VECTORS:
        raise ValueError(
            f"expected {EXPECTED_DISTINCT_VECTORS} distinct latent vectors, found {distinct}"
        )

    shared = [
        sorted(group["Key"].to_list())
        for _, group in frame.group_by(["VectorID"])
        if group.height > 1
    ]
    expected_shared = [["galactose", "mannose"], ["lactose", "maltose"]]
    if sorted(shared) != expected_shared:
        raise ValueError(f"colliding compounds {sorted(shared)} != {expected_shared}")


def main() -> None:
    """Entry point."""
    project_dir = Path(os.environ["PROJECT_DIR"])
    data_dir = Path(os.environ["DATA_DIR"]) / "bloom_256"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha", default="0.5", help="trinarization threshold directory")
    parser.add_argument(
        "--chemistry",
        type=Path,
        default=project_dir / "results/condition_chemistry.csv",
        help="canonical compound table",
    )
    parser.add_argument(
        "--out-dir", type=Path, default=project_dir / "results/latent_space"
    )
    args = parser.parse_args()

    chemistry = pl.read_csv(args.chemistry)
    mapping = label_to_key(chemistry)

    records, per_dataset_vectors, conflicts = collect(
        data_dir / f"sigma_{args.alpha}", mapping
    )

    unreached = sorted(set(chemistry["Key"].to_list()) - set(records))
    if unreached:
        raise ValueError(f"compound keys never seen in any feather: {unreached}")

    frame = build_frame(records, chemistry)
    check(frame)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    frame.write_csv(args.out_dir / "latent_conditions.csv")
    frame.with_columns(pl.col("Class").cast(pl.Categorical)).write_parquet(
        args.out_dir / "latent_conditions.parquet"
    )
    collisions = audit(frame, conflicts)
    collisions.write_csv(args.out_dir / "latent_collisions.csv")

    print(f"{frame.height} compounds, {frame['VectorID'].n_unique()} distinct latent vectors")
    print("  classes: " + ", ".join(f"{k}={v}" for k, v in sorted(EXPECTED_CLASS_SIZES.items())))
    cross_check(project_dir, per_dataset_vectors)
    print(f"  {collisions.height} collision rows")
    print(f"-> {args.out_dir}")


if __name__ == "__main__":
    main()
