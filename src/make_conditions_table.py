"""Build the growth-condition chemistry table (``results/conditions.tex``).

The three Bloom panels were assayed in overlapping but not identical sets of conditions,
and they do not spell them the same way: Bloom2013 writes ``copper`` where the other two
write ``CuSO4``, and Bloom2015 capitalises names the others leave lower case. Anything
that reports per-condition results therefore needs a single reconciled list, which is what
this script produces.

Three artefacts, all regenerated from scratch on every run:

``results/conditions.tex``
    A ``longtable`` of compound name, 2D structure, SMILES and dataset membership. Compiles
    on its own; the block between ``BEGIN TABLE`` and ``END TABLE`` can be lifted into the
    paper.
``results/condition_chemistry.csv``
    The same content without the markup -- name, PubChem CID, formula, SMILES and one
    boolean column per dataset -- for anything that needs the mapping programmatically.
``results/condition_structures/``
    One depiction per compound. PubChem PNGs by default, so the table compiles with no
    chemistry toolkit installed; ``--rdkit`` redraws them as vector PDFs instead, which
    ``\\includegraphics`` picks up in preference to the PNGs without any edit to the tex.

Membership is read from the ``Condition`` column of the feathers themselves rather than
transcribed, so the table cannot drift from the data. The chemistry is fetched once from
PubChem and cached in ``results/condition_chemistry.csv``; delete that file to refetch.

    python src/make_conditions_table.py
    python src/make_conditions_table.py --rdkit     # vector structures, needs rdkit
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import urllib.parse
import urllib.request
from pathlib import Path

import polars as pl
from dotenv import load_dotenv
from rich.console import Console

load_dotenv()

console = Console()

DATASETS = {  # display name -> feather stem
    "Bloom2013": "bloom2013",
    "Bloom2015": "bloom2015",
    "Bloom2019": "bloom2019_BYxRM",
}

# Canonical key -> the PubChem query that identifies the compound. A bare CID is used where
# the name is ambiguous (PubChem's `cisplatin` record carries no geometry, `cobalt chloride`
# does not resolve by that name, and tunicamycin is a homologue mixture with no single
# record) -- see the footnotes in the generated table.
CHEMISTRY = {
    "4NQO": ("name", "4-nitroquinoline 1-oxide", "4-Nitroquinoline 1-oxide"),
    "6-azauracil": ("name", "6-azauracil", "6-Azauracil"),
    "berbamine": ("name", "berbamine", "Berbamine"),
    "caffeine": ("name", "caffeine", "Caffeine"),
    "CaCl2": ("name", "calcium chloride", "Calcium chloride"),
    "CdCl2": ("name", "cadmium chloride", "Cadmium chloride"),
    "cisplatin": ("cid", "5460033", "Cisplatin"),
    "CoCl2": ("cid", "3032536", "Cobalt(II) chloride"),
    "congo_red": ("name", "congo red", "Congo red"),
    "CuSO4": ("name", "copper(II) sulfate", "Copper(II) sulfate"),
    "cycloheximide": ("name", "cycloheximide", "Cycloheximide"),
    "diamide": ("name", "diamide", "Diamide"),
    "egta": ("name", "egta", "EGTA"),
    "ethanol": ("name", "ethanol", "Ethanol"),
    "fluconazole": ("name", "fluconazole", "Fluconazole"),
    "fluorocytosine": ("name", "flucytosine", "Flucytosine (5-fluorocytosine)"),
    "fluorouracil": ("name", "5-fluorouracil", "5-Fluorouracil"),
    "formamide": ("name", "formamide", "Formamide"),
    "fructose": ("name", "D-fructose", "D-Fructose"),
    "galactose": ("name", "D-galactose", "D-Galactose"),
    "glycerol": ("name", "glycerol", "Glycerol"),
    "H2O2": ("name", "hydrogen peroxide", "Hydrogen peroxide"),
    "hydroquinone": ("name", "hydroquinone", "Hydroquinone"),
    "hydroxybenzaldehyde": ("name", "4-hydroxybenzaldehyde", "4-Hydroxybenzaldehyde"),
    "hydroxyurea": ("name", "hydroxyurea", "Hydroxyurea"),
    "indoleacetic_acid": ("name", "indole-3-acetic acid", "Indole-3-acetic acid"),
    "lactate": ("name", "L-lactic acid", "L-Lactic acid"),
    "lactose": ("name", "lactose", "Lactose"),
    "LiCl": ("name", "lithium chloride", "Lithium chloride"),
    "maltose": ("name", "maltose", "Maltose"),
    "mannose": ("name", "D-mannose", "D-Mannose"),
    "menadione": ("name", "menadione", "Menadione (vitamin K3)"),
    "methotrexate": ("name", "methotrexate", "Methotrexate"),
    "MgCl2": ("name", "magnesium chloride", "Magnesium chloride"),
    "MgSO4": ("name", "magnesium sulfate", "Magnesium sulfate"),
    "MnSO4": ("name", "manganese(II) sulfate", "Manganese(II) sulfate"),
    "neomycin": ("name", "neomycin B", "Neomycin B"),
    "paraquat": ("name", "paraquat dichloride", "Paraquat dichloride"),
    "raffinose": ("name", "raffinose", "Raffinose"),
    "SDS": ("name", "sodium dodecyl sulfate", "Sodium dodecyl sulfate"),
    "sorbitol": ("name", "D-sorbitol", "D-Sorbitol"),
    "sucrose": ("name", "sucrose", "Sucrose"),
    "trehalose": ("cid", "7427", "alpha,alpha-Trehalose"),
    "tunicamycin": ("cid", "11104835", "Tunicamycin B2"),
    "xylose": ("name", "D-xylose", "D-Xylose"),
    "zeocin": ("name", "zeocin", "Zeocin (phleomycin D1)"),
}

# Condition labels that name the same treatment as a differently-keyed entry above.
ALIAS = {"copper": "cuso4"}

# Footnote markers, keyed by canonical name. The text lives in TEMPLATE.
FOOTNOTES = {
    "CuSO4": "a", "cisplatin": "b", "tunicamycin": "c",
    "neomycin": "d", "zeocin": "d", "lactate": "e",
}

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
PROPERTIES = "MolecularFormula,MolecularWeight,SMILES,IUPACName,Title"


def normalise(label: str) -> str:
    """Reduce a raw condition label to a form comparable across datasets.

    Parameters
    ----------
    label : str
        A value from a dataset's ``Condition`` column.

    Returns
    -------
    str
        Lower case with separators removed, then passed through :data:`ALIAS`.
    """
    key = re.sub(r"[\s_-]", "", label).lower()
    return ALIAS.get(key, key)


def read_membership(data_dir: Path) -> dict[str, dict[str, str]]:
    """Which datasets contain each compound, and under which label.

    Parameters
    ----------
    data_dir : Path
        Directory holding ``<stem>_max.feather`` for each dataset.

    Returns
    -------
    dict[str, dict[str, str]]
        Canonical key -> dataset -> the raw label that dataset uses.

    Raises
    ------
    KeyError
        If a condition appears in the data with no entry in :data:`CHEMISTRY`. Silently
        dropping it would understate a dataset, so this stops instead.
    """
    by_norm = {normalise(key): key for key in CHEMISTRY}

    membership: dict[str, dict[str, str]] = {}
    for dataset, stem in DATASETS.items():
        labels = (
            pl.read_ipc(data_dir / f"{stem}_max.feather", columns=["Condition"])["Condition"]
            .unique()
            .to_list()
        )
        for label in labels:
            key = by_norm.get(normalise(label))
            if key is None:
                raise KeyError(
                    f"{dataset} condition {label!r} has no CHEMISTRY entry "
                    f"(normalised to {normalise(label)!r})"
                )
            membership.setdefault(key, {})[dataset] = label
        console.log(f"{dataset}: {len(labels)} conditions")

    unused = set(CHEMISTRY) - set(membership)
    if unused:
        console.log(f"[yellow]{len(unused)} CHEMISTRY entries unused:[/yellow] {sorted(unused)}")

    return membership


def fetch_chemistry(cache_path: Path) -> dict[str, dict]:
    """Structure identifiers for every compound, from the cache or from PubChem.

    Parameters
    ----------
    cache_path : Path
        ``results/condition_chemistry.csv``. Read if present; otherwise every compound is
        fetched and the file written.

    Returns
    -------
    dict[str, dict]
        Canonical key -> record with ``CID``, ``MolecularFormula``, ``SMILES``.
    """
    if cache_path.exists():
        rows = list(csv.DictReader(cache_path.open()))
        cached = {r["Key"]: r for r in rows if r["Key"] in CHEMISTRY}
        if set(cached) == set(CHEMISTRY):
            console.log(f"chemistry from cache ({len(cached)} compounds)")
            return cached
        console.log("[yellow]cache incomplete, refetching[/yellow]")

    records = {}
    for key, (kind, query, _) in CHEMISTRY.items():
        url = f"{PUBCHEM}/{kind}/{urllib.parse.quote(query)}/property/{PROPERTIES}/JSON"
        with urllib.request.urlopen(url, timeout=60) as handle:
            record = json.load(handle)["PropertyTable"]["Properties"][0]
        records[key] = {
            "Key": key,
            "CID": str(record["CID"]),
            "MolecularFormula": record["MolecularFormula"],
            "SMILES": record["SMILES"],
        }
        console.log(f"{key:22s} CID {record['CID']:>10}  {record['MolecularFormula']}")

    return records


def fetch_depictions(records: dict[str, dict], out_dir: Path) -> None:
    """Download one PubChem 2D depiction per compound, skipping any already present."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for key, record in sorted(records.items()):
        dest = out_dir / f"{key}.png"
        if dest.exists() and dest.stat().st_size > 500:
            continue
        url = f"{PUBCHEM}/cid/{record['CID']}/PNG?image_size=500x500"
        with urllib.request.urlopen(url, timeout=60) as handle:
            dest.write_bytes(handle.read())
        console.log(f"depiction {key}")


def draw_with_rdkit(records: dict[str, dict], out_dir: Path) -> None:
    """Redraw every structure as a vector PDF from its SMILES.

    ``\\includegraphics`` is given no extension in the table, so a ``.pdf`` written here is
    picked up ahead of the ``.png`` that :func:`fetch_depictions` downloaded.

    Parameters
    ----------
    records : dict[str, dict]
        Output of :func:`fetch_chemistry`.
    out_dir : Path
        ``results/condition_structures``.
    """
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D

    out_dir.mkdir(parents=True, exist_ok=True)
    for key, record in sorted(records.items()):
        mol = Chem.MolFromSmiles(record["SMILES"])
        if mol is None:
            console.log(f"[red]{key}: rdkit could not parse the SMILES[/red]")
            continue
        Draw.rdDepictor.Compute2DCoords(mol)
        canvas = rdMolDraw2D.MolDraw2DCairo(500, 500)
        canvas.drawOptions().clearBackground = False
        rdMolDraw2D.PrepareAndDrawMolecule(canvas, mol)
        canvas.FinishDrawing()
        (out_dir / f"{key}.pdf").write_bytes(canvas.GetDrawingText())
        console.log(f"drew {key}")


def escape(text: str) -> str:
    """Escape the characters LaTeX treats specially."""
    table = {
        "\\": r"\textbackslash{}", "#": r"\#", "%": r"\%", "&": r"\&", "_": r"\_",
        "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
        "$": r"\$",
    }
    return "".join(table.get(char, char) for char in text)


def wrap_smiles(smiles: str, chunk: int = 6) -> str:
    """Escape a SMILES and permit a line break every ``chunk`` characters.

    A SMILES is one unbroken word, so TeX has nowhere to break it and it runs off the
    column. ``\\allowbreak`` inserts a zero-width, hyphen-free break opportunity.
    """
    pieces = [escape(char) for char in smiles]
    out = []
    for index, piece in enumerate(pieces, start=1):
        out.append(piece)
        if index % chunk == 0:
            out.append(r"\allowbreak{}")
    return "".join(out)


def sort_key(display: str) -> str:
    """Alphabetise ignoring locants and stereo-descriptor prefixes.

    ``4-Hydroxybenzaldehyde`` files under H and ``alpha,alpha-Trehalose`` under T, which is
    how a reader scanning the table will look for them.
    """
    name = display
    while True:
        stripped = re.sub(r"^(\d+|[DLNRS]|alpha,alpha|alpha|beta|cis|trans)-", "", name)
        if stripped == name:
            return name.lower()
        name = stripped


def latex_name(display: str) -> str:
    """Typeset a display name, italicising stereo descriptors and setting subscripts."""
    name = escape(display)
    name = name.replace("alpha,alpha-", r"$\alpha,\alpha$-")
    name = re.sub(r"^([DL])-", r"\\textsc{\1}-", name)
    return name.replace("K3", "K$_3$")


def build_tex(membership: dict[str, dict[str, str]], records: dict[str, dict]) -> str:
    """Render the whole ``conditions.tex`` document."""
    order = sorted(membership, key=lambda key: sort_key(CHEMISTRY[key][2]))

    rows = []
    for key in order:
        name = latex_name(CHEMISTRY[key][2])
        if key in FOOTNOTES:
            name += r"\textsuperscript{" + FOOTNOTES[key] + "}"
        labels = []
        for dataset in DATASETS:
            label = membership[key].get(dataset)
            if label and label not in labels:
                labels.append(label)
        labels_tex = ", ".join(r"\texttt{" + escape(label) + "}" for label in labels)
        ticks = " & ".join(
            r"\yes" if dataset in membership[key] else r"\no" for dataset in DATASETS
        )
        rows.append(
            f"{name}\\newline {{\\scriptsize {labels_tex}}}\n"
            f"  & \\structure{{{key}}}\n"
            f"  & \\smi{{{wrap_smiles(records[key]['SMILES'])}}}\n"
            f"  & {ticks} \\\\"
        )

    counts = {ds: sum(ds in m for m in membership.values()) for ds in DATASETS}
    shared = sum(len(m) == len(DATASETS) for m in membership.values())

    return (
        TEMPLATE.replace("__BODY__", "\n".join(rows))
        .replace("__TOTAL__", str(len(order)))
        .replace("__SHARED__", str(shared))
        .replace("__N2013__", str(counts["Bloom2013"]))
        .replace("__N2015__", str(counts["Bloom2015"]))
        .replace("__N2019__", str(counts["Bloom2019"]))
    )


def write_csv(
    membership: dict[str, dict[str, str]], records: dict[str, dict], path: Path
) -> None:
    """Write the markup-free version of the same table."""
    fields = ["Key", "Name", "CID", "MolecularFormula", "SMILES"]
    fields += list(DATASETS) + [f"Label{ds}" for ds in DATASETS]

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for key in sorted(membership, key=lambda k: sort_key(CHEMISTRY[k][2])):
            row = {
                "Key": key,
                "Name": CHEMISTRY[key][2],
                "CID": records[key]["CID"],
                "MolecularFormula": records[key]["MolecularFormula"],
                "SMILES": records[key]["SMILES"],
            }
            for dataset in DATASETS:
                row[dataset] = int(dataset in membership[key])
                row[f"Label{dataset}"] = membership[key].get(dataset, "")
            writer.writerow(row)


TEMPLATE = r"""%% results/conditions.tex -- the growth conditions assayed in the three Bloom panels.
%%
%% Generated by src/make_conditions_table.py; do not hand-edit. Compound membership comes
%% from the Condition column of bloom{2013,2015,2019_BYxRM}_max.feather, so the table
%% cannot drift from the data; structures and SMILES come from PubChem, and the CIDs are
%% recorded in results/condition_chemistry.csv.
%%
%% Compiles on its own (pdflatex conditions.tex) and reads its images from
%% results/condition_structures/. To use it in the paper, copy the block between BEGIN
%% TABLE and END TABLE and load the packages named just below BEGIN TABLE.

\documentclass[10pt,a4paper]{article}
\usepackage[margin=1.8cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{graphicx}
\usepackage{array}
\usepackage{longtable}
\usepackage{booktabs}
\usepackage{amssymb}
\usepackage{xcolor}

\begin{document}

%% ------------------------------------------------------------------ BEGIN TABLE
%% requires: graphicx, array, longtable, booktabs, amssymb, xcolor
%%
%% \structure takes no file extension, so a vector .pdf written by
%% `python src/make_conditions_table.py --rdkit` is used in preference to the .png
%% without any change here.
\providecommand{\structure}[1]{%
  \includegraphics[width=2.3cm]{condition_structures/#1}}
\providecommand{\smi}[1]{{\ttfamily\fontsize{5}{6}\selectfont #1}}
\providecommand{\yes}{$\checkmark$}
\providecommand{\no}{\textcolor{black!30}{--}}

\begingroup
\setlength{\tabcolsep}{4pt}
\renewcommand{\arraystretch}{1.15}
\footnotesize
\begin{longtable}{@{}>{\raggedright\arraybackslash}m{3.2cm}
                    >{\centering\arraybackslash}m{2.5cm}
                    >{\raggedright\arraybackslash}m{7.2cm}
                    c c c@{}}
\caption{The __TOTAL__ growth conditions assayed across the three yeast segregant panels,
__SHARED__ of them common to all three. The typewriter strings beneath each compound name
are the exact values taken by the \texttt{Condition} column of that dataset, which differ
between panels in case and in naming.\label{tab:conditions}}\\
\toprule
Compound & Structure & SMILES
  & \rotatebox{90}{Bloom2013\,} & \rotatebox{90}{Bloom2015\,}
  & \rotatebox{90}{Bloom2019\textsuperscript{f}\,} \\
\midrule
\endfirsthead
\caption[]{\emph{(continued)}}\\
\toprule
Compound & Structure & SMILES
  & \rotatebox{90}{Bloom2013\,} & \rotatebox{90}{Bloom2015\,}
  & \rotatebox{90}{Bloom2019\textsuperscript{f}\,} \\
\midrule
\endhead
\midrule
\multicolumn{6}{r@{}}{\emph{continued on next page}}\\
\endfoot
\midrule
\multicolumn{3}{@{}l}{\textbf{Conditions per dataset}} & __N2013__ & __N2015__ & __N2019__ \\
\bottomrule
\endlastfoot
__BODY__
\end{longtable}
\endgroup

\noindent\footnotesize
\textsuperscript{a}~Bloom2013 labels this condition \texttt{copper} and the other two
panels label the same treatment \texttt{CuSO4}; it is counted once.\\
\textsuperscript{b}~SMILES does not encode square-planar geometry, so the string shown
does not distinguish cisplatin from its clinically inactive \emph{trans} isomer.\\
\textsuperscript{c}~Tunicamycin is supplied as a mixture of homologues differing in the
length and branching of the fatty-acyl chain; the B2 homologue is drawn as
representative.\\
\textsuperscript{d}~Neomycin and zeocin are likewise mixtures of related congeners; the
principal component is drawn (neomycin B, phleomycin D1).\\
\textsuperscript{e}~Used as a carbon source in the form of its sodium salt; the free acid
is drawn.\\
\textsuperscript{f}~Bloom2019 is the BYxRM cross, the only 2019 panel carried through this
analysis.
%% ------------------------------------------------------------------ END TABLE

\end{document}
"""


def main() -> None:
    """Entry point."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--rdkit",
        action="store_true",
        help="redraw the structures as vector PDFs from the SMILES (needs rdkit)",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(os.environ["DATA_DIR"]) / "bloom_256" / "sigma_0.5",
        help="directory holding the *_max.feather files",
    )
    args = parser.parse_args()

    project = Path(os.environ["PROJECT_DIR"])
    results = project / "results"
    csv_path = results / "condition_chemistry.csv"

    membership = read_membership(args.data_dir)
    records = fetch_chemistry(csv_path)

    write_csv(membership, records, csv_path)
    fetch_depictions(records, results / "condition_structures")
    if args.rdkit:
        draw_with_rdkit(records, results / "condition_structures")

    tex_path = results / "conditions.tex"
    tex_path.write_text(build_tex(membership, records))
    console.log(f"[green]done[/green] -> {tex_path} ({len(membership)} compounds)")


if __name__ == "__main__":
    main()
