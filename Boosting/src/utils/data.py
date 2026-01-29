# A utility file containing functions to print statistics about the data

import polars as pl
import polars.selectors as cs
from pathlib import Path
from typing import Union

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils import get_data


def get_data_stats(data_path: Union[Path, str], run_type: str = "full") -> None:
    """Prints statistics about the data.

    Parameters
    ----------
    data_path : Union[Path, str]
        Path to the data file.
    run_type: str
        One of `full`, `geno_only`, `chem_only`
    """
    df = get_data(data_path=data_path, run_type=run_type)
    print(f"Shape of the data: {df.shape}")
    print(f"Number of unique Strains: \n {df['Strain'].n_unique()}")
    print(f"Number of unique Conditions: \n {df['Condition'].value_counts()}")
    print(f"Number of unique Phenotypes: \n {df['Phenotype'].value_counts()}")

    # A little more involved: identifying how many genes have some variation

    gene_df = df.select("Strain", cs.starts_with("Y")).unique()
    gene_df = gene_df.melt(
        id_vars="Strain", variable_name="Gene", value_name="Mutation"
    )
    gene_df = gene_df.group_by("Gene").agg(pl.col("Mutation").var().alias("Variation"))

    print(gene_df.filter(pl.col("Variation") > 0))


if __name__ == "__main__":
    # Small application to get some statistics about the data

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "--data_path",
        type=str,
    )
    parser.add_argument("--run_type", type=str, default="full")
    args = parser.parse_args()

    get_data_stats(data_path=args.data_path, run_type=args.run_type)
