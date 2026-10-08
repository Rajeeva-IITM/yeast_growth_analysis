"""Build the per-condition supplementary table (``results/condition_summary/``).

``results/conditions.tex`` says *which* compounds were assayed and what they are
chemically. It does not say how much data each one contributes, which is what a reader
needs in order to weigh a per-condition result. This script adds that layer: for every
(dataset, condition) pair it reports the study the measurements come from, how many
segregants were assayed, how many survived trinarisation, and how the survivors split
between the high- and low-growth classes.

Three numbers are kept deliberately distinct, because they are not interchangeable:

``SegregantsAssayed``
    Segregants with a non-missing growth measurement for that condition in the published
    phenotype file, restricted to the segregants that made it into the modelling data.
``MeasurementsAssayed``
    Non-missing measurements of that condition. Equal to ``SegregantsAssayed`` everywhere
    except Bloom2015, whose phenotype release reports two replicate columns per condition
    and so measures some segregants twice.
``SegregantsUsed``
    Segregants that contribute at least one row to the modelling data for that condition,
    i.e. whose measurement fell outside the +/- alpha*sd dead zone.
``Observations``
    Rows in the modelling data, equal to ``LowGrowth + HighGrowth``. This is the number
    directly comparable to ``MeasurementsAssayed``; comparing it to a segregant count is
    an error for Bloom2015.
``ReplicateAgreement``
    Where a segregant was measured twice and both measurements survived trinarisation,
    the fraction of such pairs that landed in the same class. Bloom2015 only.

``Phenotype`` in the source parquets is trinarised at mean +/- alpha*sd with missing cells
mean-imputed, so imputed measurements always land in the discarded middle class. The
modelling feathers keep only the two extreme classes, coded 0 = low growth, 1 = high
growth. The gap between ``SegregantsAssayed`` and ``SegregantsUsed`` is therefore the
middle class, and is expected to be large -- roughly 40% of measurements at alpha 0.5.

Outputs, regenerated from scratch on every run:

``results/condition_summary/condition_summary.csv``
    One row per (dataset, condition); 88 rows at the paper's three panels.
``results/condition_summary/condition_summary.tex``
    The same content as a standalone ``longtable`` supplement, grouped by compound.

Nothing is transcribed: membership, counts and class balance are all read from the files
the models were trained on, and the label-to-phenotype-column mapping is rebuilt from the
dataset builder's own ``name_map`` configs.

    pixi run python src/make_condition_summary.py
    pixi run python src/make_condition_summary.py --alpha 1.0
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
from pathlib import Path

import polars as pl
from dotenv import load_dotenv
from rich.console import Console

from make_conditions_table import CHEMISTRY, escape, latex_name, normalise, sort_key
from make_latent_dataset import CLASSES

load_dotenv()

console = Console()

# Display name, feather stem, phenotype file and builder config for each panel. `id_column`
# differs because the Bloom2019 phenotype release keys on `id` and holds all three crosses
# in one file; the segregant names are disjoint between crosses, so selecting the rows of
# one cross by name is unambiguous (verified: 0 shared names).
DATASETS = {
    "Bloom2013": {
        "stem": "bloom2013",
        "source": "Bloom et al. 2013",
        "cross": "BY x RM",
        "phenotype": "bloom2013_pheno.tsv",
        "separator": "\t",
        "id_column": "Strain",
        "config": "bloom2013.json",
    },
    "Bloom2015": {
        "stem": "bloom2015",
        "source": "Bloom et al. 2015",
        "cross": "BY x RM",
        "phenotype": "bloom2015_pheno.tsv",
        "separator": ",",
        "id_column": "Strain",
        "config": "bloom2015.json",
    },
    "Bloom2019_BYxRM": {
        "stem": "bloom2019_BYxRM",
        "source": "Bloom et al. 2019",
        "cross": "BY x RM",
        "phenotype": "bloom2019_pheno.tsv",
        "separator": "\t",
        "id_column": "id",
        "config": "BYxRM.json",
    },
}

CLASS_LABEL = {
    "salt": "Salt",
    "carbon_source": "Carbon source",
    "other": "Other",
}

# Media columns present in the phenotype releases that carry no compound and are excluded
# from the datasets by the builder. Listed so the table can say so rather than leave the
# reader to notice the arithmetic does not close.
CONTROL_COLUMNS = {
    "Bloom2015": ["YNB", "YPD"],
    "Bloom2019_BYxRM": ["YNB", "YPD", "EtOH_Glucose"],
}

EXPECTED_CONDITIONS = {"Bloom2013": 39, "Bloom2015": 18, "Bloom2019_BYxRM": 31}
EXPECTED_COMPOUNDS = 46


def class_of() -> dict[str, str]:
    """Canonical key -> chemical class, inverted from :data:`CLASSES`."""
    lookup = {key: name for name, members in CLASSES.items() for key in members}
    missing = set(CHEMISTRY) - set(lookup)
    if missing:
        raise KeyError(f"compounds with no class: {sorted(missing)}")
    return lookup


def phenotype_columns(
    dataset: str, header: list[str], name_map: dict[str, str]
) -> dict[str, dict]:
    """Map each canonical compound key to its columns in the phenotype file.

    Parameters
    ----------
    dataset : str
        Key of :data:`DATASETS`.
    header : list[str]
        Column names of the phenotype file, first (identifier) column removed.
    name_map : dict[str, str]
        The builder's ``name_map``, applied before normalisation so that e.g. Bloom2015's
        ``CopperSulfate`` and Bloom2019's ``Copper_Sulfate`` both reach ``CuSO4``.

    Returns
    -------
    dict[str, dict]
        Canonical key -> ``{"columns": [...], "dose": str}``. Control media are dropped.

    Notes
    -----
    Bloom2015 names its columns ``<condition> _replicate_<n>``; Bloom2019 names them
    ``<condition>;<dose>;<batch>``. Both are reduced to the bare condition here, and the
    dose is kept where the file reports one.
    """
    by_norm = {normalise(key): key for key in CHEMISTRY}

    grouped: dict[str, dict] = {}
    for column in header:
        base, dose = column, ""
        if dataset == "Bloom2015":
            base = re.sub(r"\s*_?replicate_\d+$", "", column).strip()
        elif dataset.startswith("Bloom2019"):
            parts = column.split(";")
            base, dose = parts[0], parts[1] if len(parts) > 1 else ""

        key = by_norm.get(normalise(name_map.get(base, base)))
        if key is None:
            if base not in CONTROL_COLUMNS.get(dataset, []):
                raise KeyError(f"{dataset} phenotype column {column!r} is unaccounted for")
            continue

        record = grouped.setdefault(key, {"columns": [], "dose": dose})
        record["columns"].append(column)
        if dose and not record["dose"]:
            record["dose"] = dose

    return grouped


def assayed(frame: pl.DataFrame, columns: list[str]) -> tuple[int, int]:
    """Count non-missing measurements of one condition.

    Parameters
    ----------
    frame : pl.DataFrame
        The phenotype release, already restricted to the segregants of one dataset.
    columns : list[str]
        Every column of ``frame`` reporting that condition -- more than one where the
        release carries replicates.

    Returns
    -------
    tuple[int, int]
        Segregants with at least one measurement, and measurements in total.
    """
    flags = [pl.col(column).is_not_null() for column in columns]
    segregants = int(frame.select(pl.any_horizontal(flags).sum()).item())
    measurements = int(frame.select(pl.sum_horizontal(flags).sum()).item())
    return segregants, measurements


def replicate_agreement(frame: pl.DataFrame) -> dict[str, float]:
    """Per condition, how often a twice-measured segregant lands in the same class.

    Parameters
    ----------
    frame : pl.DataFrame
        ``Strain``, ``Condition`` and ``Phenotype`` for one dataset.

    Returns
    -------
    dict[str, float]
        Condition -> agreeing fraction, over the segregants contributing two rows to that
        condition. Empty for datasets without replicates.

    Notes
    -----
    Only pairs in which *both* replicates survived trinarisation are counted, so this is
    the concordance among the measurements a model actually sees, not the reproducibility
    of the underlying growth assay.
    """
    pairs = (
        frame.group_by(["Strain", "Condition"])
        .agg(pl.len().alias("n"), pl.col("Phenotype").n_unique().alias("classes"))
        .filter(pl.col("n") == 2)
    )
    if pairs.height == 0:
        return {}

    return {
        row["Condition"]: row["Agreement"]
        for row in pairs.group_by("Condition")
        .agg(((pl.col("classes") == 1).sum() / pl.len()).alias("Agreement"))
        .to_dicts()
    }


def summarise(data_dir: Path, pheno_dir: Path, config_dir: Path) -> list[dict]:
    """Assemble one record per (dataset, condition).

    Parameters
    ----------
    data_dir : Path
        Directory holding ``<stem>_max.feather`` for the chosen threshold.
    pheno_dir : Path
        Directory holding the published phenotype tables.
    config_dir : Path
        Directory holding the dataset builder's json configs.

    Returns
    -------
    list[dict]
        Records sorted by compound, then by dataset in :data:`DATASETS` order.
    """
    classes = class_of()
    records = []

    for dataset, spec in DATASETS.items():
        frame = pl.read_ipc(
            data_dir / f"{spec['stem']}_max.feather",
            columns=["Strain", "Condition", "Phenotype"],
        )
        strains = frame["Strain"].unique()

        used = (
            frame.group_by("Condition")
            .agg(
                pl.len().alias("Observations"),
                pl.col("Strain").n_unique().alias("SegregantsUsed"),
                (pl.col("Phenotype") == 1).sum().alias("HighGrowth"),
                (pl.col("Phenotype") == 0).sum().alias("LowGrowth"),
            )
            .to_dicts()
        )
        used = {row.pop("Condition"): row for row in used}
        agreement = replicate_agreement(frame)

        if len(used) != EXPECTED_CONDITIONS[dataset]:
            raise ValueError(
                f"{dataset}: {len(used)} conditions, expected {EXPECTED_CONDITIONS[dataset]}"
            )

        pheno = pl.read_csv(
            pheno_dir / spec["phenotype"],
            separator=spec["separator"],
            null_values=["NA", ""],
            infer_schema_length=None,
        ).filter(pl.col(spec["id_column"]).is_in(strains))

        if pheno.height != len(strains):
            raise ValueError(
                f"{dataset}: {pheno.height} of {len(strains)} segregants found in "
                f"{spec['phenotype']}"
            )

        name_map = json.loads((config_dir / spec["config"]).read_text()).get("name_map", {})
        columns = phenotype_columns(
            dataset, [c for c in pheno.columns if c != spec["id_column"]], name_map
        )

        by_norm = {normalise(key): key for key in CHEMISTRY}
        for label, counts in used.items():
            key = by_norm[normalise(label)]
            if key not in columns:
                raise KeyError(f"{dataset} condition {label!r} has no phenotype column")

            total = counts["HighGrowth"] + counts["LowGrowth"]
            n_segregants, n_measurements = assayed(pheno, columns[key]["columns"])
            records.append(
                {
                    "Compound": CHEMISTRY[key][2],
                    "Key": key,
                    "Class": CLASS_LABEL[classes[key]],
                    "Dataset": dataset,
                    "Source": spec["source"],
                    "Cross": spec["cross"],
                    "Label": label,
                    "Dose": columns[key]["dose"],
                    "Replicates": len(columns[key]["columns"]),
                    "SegregantsAssayed": n_segregants,
                    "MeasurementsAssayed": n_measurements,
                    "SegregantsUsed": counts["SegregantsUsed"],
                    "Observations": counts["Observations"],
                    "LowGrowth": counts["LowGrowth"],
                    "HighGrowth": counts["HighGrowth"],
                    "PercentHigh": round(100 * counts["HighGrowth"] / total, 1),
                    "PercentRetained": round(100 * total / n_measurements, 1),
                    "ReplicateAgreement": (
                        round(agreement[label], 3) if label in agreement else ""
                    ),
                }
            )

        console.log(
            f"{dataset}: {len(used)} conditions, {len(strains)} segregants, "
            f"{frame.height} rows"
        )

    compounds = {record["Key"] for record in records}
    if len(compounds) != EXPECTED_COMPOUNDS:
        raise ValueError(f"{len(compounds)} compounds, expected {EXPECTED_COMPOUNDS}")

    order = list(DATASETS)
    records.sort(key=lambda r: (sort_key(r["Compound"]), order.index(r["Dataset"])))
    return records


FIELDS = [
    "Compound", "Key", "Class", "Dataset", "Source", "Cross", "Label", "Dose",
    "Replicates", "SegregantsAssayed", "MeasurementsAssayed", "SegregantsUsed",
    "Observations", "LowGrowth", "HighGrowth", "PercentHigh", "PercentRetained",
    "ReplicateAgreement",
]


def write_csv(records: list[dict], path: Path) -> None:
    """Write the machine-readable version of the table."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(records)


def dataset_totals(records: list[dict]) -> dict[str, dict]:
    """Column sums per dataset, for the summary block at the foot of the table."""
    totals = {}
    for dataset in DATASETS:
        rows = [r for r in records if r["Dataset"] == dataset]
        observations = sum(r["Observations"] for r in rows)
        high = sum(r["HighGrowth"] for r in rows)
        measurements = sum(r["MeasurementsAssayed"] for r in rows)
        totals[dataset] = {
            "Conditions": len(rows),
            "SegregantsAssayed": sum(r["SegregantsAssayed"] for r in rows),
            "MeasurementsAssayed": measurements,
            "SegregantsUsed": sum(r["SegregantsUsed"] for r in rows),
            "Observations": observations,
            "LowGrowth": observations - high,
            "HighGrowth": high,
            "PercentHigh": round(100 * high / observations, 1),
            "PercentRetained": round(100 * observations / measurements, 1),
        }
    return totals


def format_dose(dose: str) -> str:
    """Typeset a dose string: thin space before the unit, micro sign where meant.

    Parameters
    ----------
    dose : str
        As recorded in the Bloom2019 phenotype header, e.g. ``75uM`` or ``5mg/mL``.

    Returns
    -------
    str
        LaTeX, or an en dash where the source release records no dose.

    Notes
    -----
    The recorded doses use only digits, letters, ``.``, ``/`` and ``%``, so ``%`` is the
    only character needing an escape and the substitutions below cannot collide with one.
    """
    if not dose:
        return r"\textendash"
    text = dose.replace("%", r"\%")
    text = re.sub(r"^([\d.]+)", r"\1\\,", text)
    return re.sub(r"\bu(?=[MgmL])", r"$\\mu$", text)


def build_tex(records: list[dict], alpha: float) -> str:
    """Render the standalone supplementary longtable."""
    order = list(DATASETS)
    body = []
    previous = None
    for record in records:
        first = record["Compound"] != previous
        if first and previous is not None:
            body.append(r"\addlinespace[2pt]")
        name = latex_name(record["Compound"]) if first else ""
        klass = record["Class"] if first else ""
        dose = format_dose(record["Dose"])
        body.append(
            f"{name} & {klass} & {record['Source']} & {dose} & "
            f"{record['SegregantsAssayed']} & {record['MeasurementsAssayed']} & "
            f"{record['LowGrowth']} & {record['HighGrowth']} & "
            f"{record['PercentHigh']:.1f} \\\\"
        )
        previous = record["Compound"]

    totals = dataset_totals(records)
    summary = []
    for dataset in order:
        row = totals[dataset]
        summary.append(
            f"{DATASETS[dataset]['source']} & {row['Conditions']} & "
            f"{row['SegregantsAssayed']} & {row['MeasurementsAssayed']} & "
            f"{row['Observations']} & {row['PercentRetained']:.1f} & "
            f"{row['LowGrowth']} & {row['HighGrowth']} & "
            f"{row['PercentHigh']:.1f} \\\\"
        )

    return (
        TEMPLATE.replace("__BODY__", "\n".join(body))
        .replace("__SUMMARY__", "\n".join(summary))
        .replace("__ALPHA__", f"{alpha:.1f}")
        .replace("__PAIRS__", str(len(records)))
        .replace("__COMPOUNDS__", str(len({r["Key"] for r in records})))
        .replace("__AGREEMENT__", f"{100 * mean_agreement(records):.1f}")
    )


def mean_agreement(records: list[dict]) -> float:
    """Replicate concordance pooled over every condition that reports one."""
    values = [r["ReplicateAgreement"] for r in records if r["ReplicateAgreement"] != ""]
    return sum(values) / len(values) if values else 0.0


TEMPLATE = r"""%% results/condition_summary/condition_summary.tex -- per-condition data summary.
%%
%% Generated by src/make_condition_summary.py; do not hand-edit. Every count is read from
%% the files the models were trained on (bloom_256/sigma___ALPHA__/*_max.feather) and from
%% the published phenotype releases, so the table cannot drift from the data.
%%
%% Companion to results/conditions.tex, which gives the chemistry of the same compounds.
%% Compiles on its own; the block between BEGIN TABLE and END TABLE lifts into the paper.

\documentclass[10pt,a4paper]{article}
\usepackage[margin=1.6cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{longtable}
\usepackage{booktabs}
\usepackage{array}

\begin{document}

%% ------------------------------------------------------------------ BEGIN TABLE
%% requires: longtable, booktabs, array
\begingroup
\setlength{\tabcolsep}{4pt}
\footnotesize

\begin{longtable}{@{}llll rr rrr@{}}
\caption{\textbf{The __COMPOUNDS__ growth conditions and the data each contributes.}
Every condition assayed in the three segregant panels, with the study it comes from, the
amount of data it contributes and the resulting class balance. Growth was trinarised
within each condition at the panel mean $\pm$ __ALPHA__ standard deviations; the middle
class is discarded and the two extremes are modelled as a binary label, so
\emph{low}~$+$~\emph{high} is smaller than \emph{measurements}. A compound assayed in more
than one panel occupies one row per panel, giving __PAIRS__ panel--condition pairs in
total. \emph{Segregants} and \emph{measurements} differ only for Bloom2015, whose
phenotype release reports two replicate measurements per segregant.}
\label{tab:condition-summary} \\
\toprule
& & & & \multicolumn{2}{c}{Assayed} & \multicolumn{3}{c}{Modelled} \\
\cmidrule(lr){5-6} \cmidrule(lr){7-9}
Compound & Class & Source & Dose & Segregants & Meas. & Low & High & \% high \\
\midrule
\endfirsthead

\multicolumn{9}{@{}l}{\textit{Table \ref{tab:condition-summary}, continued}} \\
\toprule
& & & & \multicolumn{2}{c}{Assayed} & \multicolumn{3}{c}{Modelled} \\
\cmidrule(lr){5-6} \cmidrule(lr){7-9}
Compound & Class & Source & Dose & Segregants & Meas. & Low & High & \% high \\
\midrule
\endhead

\midrule
\multicolumn{9}{r@{}}{\textit{continued on the next page}} \\
\endfoot

\bottomrule
\endlastfoot

__BODY__
\end{longtable}

\begin{table}[!ht]
\centering
\footnotesize
\caption{\textbf{Per-panel totals.} Column sums over
Table~\ref{tab:condition-summary}; the assayed columns count segregant--condition pairs,
not distinct segregants. \emph{Used} is the number of rows entering the models and
\emph{kept} the percentage of measurements that survived trinarisation.}
\label{tab:condition-summary-totals}
\begin{tabular}{@{}l r rrr r rr r@{}}
\toprule
& & \multicolumn{2}{c}{Assayed} & & & \multicolumn{3}{c}{Modelled} \\
\cmidrule(lr){3-4} \cmidrule(lr){7-9}
Source & Conditions & Segregants & Meas. & Used & \% kept & Low & High & \% high \\
\midrule
__SUMMARY__
\bottomrule
\end{tabular}
\end{table}

\endgroup

%% Footnotes worth carrying into the paper:
%%  - Dose is reported only where the source phenotype release records it; the Bloom2013
%%    and Bloom2015 releases do not, and the concentrations are given in those papers.
%%  - Control media (YPD, YNB, and Bloom2019's ethanol+glucose) appear in the phenotype
%%    releases but carry no compound and are excluded from every dataset.
%%  - Missing measurements are mean-imputed before trinarisation and therefore always land
%%    in the discarded middle class, so no imputed value reaches a model.
%%  - Bloom2015 replicate pairs in which both measurements survived trinarisation agree on
%%    the class __AGREEMENT__% of the time (per-condition values in
%%    condition_summary.csv), which bounds how much of that panel's signal is measurement
%%    noise.
%% ------------------------------------------------------------------ END TABLE

\end{document}
"""


def parse_args() -> argparse.Namespace:
    """Command-line interface."""
    root = Path(os.environ.get("PROJECT_DIR", Path.cwd()))
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--alpha", type=float, default=0.5,
                        help="trinarisation threshold whose datasets are summarised")
    parser.add_argument("--data-dir", type=Path,
                        default=Path(os.environ.get("DATA_DIR", "")) / "bloom_256",
                        help="directory holding the sigma_<alpha> modelling feathers")
    parser.add_argument("--pheno-dir", type=Path,
                        default=Path(os.environ.get("DATA_DIR", "")) / "pheno",
                        help="directory holding the published phenotype releases")
    parser.add_argument("--config-dir", type=Path,
                        default=Path.home() / "projects/bloom_data/configs",
                        help="directory holding the dataset builder configs")
    parser.add_argument("--out-dir", type=Path, default=root / "results/condition_summary",
                        help="where the csv and tex are written")
    return parser.parse_args()


def main() -> None:
    """Build both artefacts and report the headline counts."""
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    records = summarise(
        args.data_dir / f"sigma_{args.alpha:.1f}", args.pheno_dir, args.config_dir
    )

    write_csv(records, args.out_dir / "condition_summary.csv")
    (args.out_dir / "condition_summary.tex").write_text(build_tex(records, args.alpha))

    totals = dataset_totals(records)
    for dataset, row in totals.items():
        console.log(
            f"{dataset:18s} {row['Conditions']:3d} conditions  "
            f"{row['Observations']:6d}/{row['MeasurementsAssayed']:6d} measurements kept "
            f"({row['PercentRetained']:.1f}%)  {row['PercentHigh']:.1f}% high"
        )
    console.log(
        f"{len(records)} panel-condition pairs over "
        f"{len({r['Key'] for r in records})} compounds -> {args.out_dir}"
    )


if __name__ == "__main__":
    main()
