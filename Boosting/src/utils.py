# utility functions

import lightgbm as lgb
import polars as pl
import polars.selectors as cs
import pickle
from pathlib import Path
from typing import Dict, List, Union, Tuple, Any
from numpy import ndarray

from scipy.stats import pearsonr, spearmanr


# Getting stuff
def get_data(
    data_path: Union[Path, str], run_type: str, return_as_Xy: bool = False
) -> Union[pl.DataFrame, Tuple[ndarray[Any, Any], ndarray[Any, Any]]]:
    """Loads the data from a given path and returns it as a DataFrame.

    Parameters
    ----------
    data_path : Union[Path, str]
        Path to the data file.

    run_type: str
        One of `full`, `geno_only`, `chem_only`

    return_as_Xy: bool
        If you want the data to be returned as [Features, Observations]

    Returns
    -------
    Union[pl.DataFrame, Tuple[pl.DataFrame, pl.DataFrame]]
        The loaded data as a DataFrame.
    """

    if isinstance(data_path, Path):
        format = data_path.suffix
    elif isinstance(data_path, str):
        format = data_path.split(".")[-1]
    else:
        raise ValueError(
            f"path must be either Path or str but it is of type {type(data_path)}"
        )

    match format:
        case "feather":
            df = pl.read_ipc(data_path)
        case "parquet":
            df = pl.read_parquet(data_path)
        case "csv":
            df = pl.read_csv(data_path)
        case _:
            raise NotImplementedError(
                "File format not supported yet. Most be one of 'feather', 'parquet', 'csv'"
            )

    if not return_as_Xy:
        # y = df['Phenotype'].to_numpy()

        match run_type:
            case "geno_only":
                df = df.select(
                    pl.col("Strain"),
                    cs.starts_with(
                        "Y"
                    ),  # Only the genotype columns, based on yeast systemic names # TODO: Hardcoded for now. But should be flexible
                    pl.col("Condition"),
                    pl.col("Phenotype"),
                )
                return df

            case "chem_only":
                df = df.select(
                    pl.col("Strain"),
                    cs.contains(
                        "latent"
                    ),  # Only the chemical information # TODO: Hardcoded for now. But should be flexible
                    pl.col("Condition"),
                    pl.col("Phenotype"),
                )
                return df
            case "full":
                return df
            case _:
                raise ValueError(
                    "Invalid run_type. Must be one of ['geno_only', 'chem_only', 'full']"
                )

    else:
        y = df["Phenotype"].to_numpy()

        match run_type:
            case "geno_only":
                X = df.select(cs.starts_with("Y")).to_numpy()
                return X, y
            case "chem_only":
                X = df.select(cs.contains("latent")).to_numpy()
                return X, y
            case "full":
                X = df.select(cs.starts_with("Y"), cs.contains("latent")).to_numpy()
                return X, y
            case _:
                raise ValueError(
                    "Invalid run_type. Must be one of ['geno_only', 'chem_only', 'full']"
                )


def get_model(model_path: Path) -> lgb.Booster:
    """Loads a LightGBM model from a pickle file and returns it.

    Parameters
    ----------
    model_path : Path
        Path to the pickle file containing the model.

    Returns
    -------
    lgb.Booster
        The loaded model.
    """
    model = pickle.load(open(model_path, "rb"))
    assert isinstance(
        model, lgb.Booster
    ), "Model is not an instance of lgb.Booster, villain."
    return model


def get_model_paths(
    run_path: str, model_type: str, suffix: str, prefix: str, model_names: List[str]
) -> Dict[str, List[Path]]:
    """Retrieves the paths of the models in a given directory based on the provided model type and
    suffix.

    Parameters:
        run_path (str): The directory containing the models.
        model_type (str): The type of model to retrieve (e.g. lgbm, dummy).
        suffix (str): The suffix of the model (e.g. full, geno_only, chem_only, dummy).
        prefix (str): The prefix of the model paths (e.g. Bloom2013_).
        model_names (List[str]): The list of models to retrieve.

    Returns:
        Dict[str, List[Path]]: A dictionary with the model names as keys and the
            paths to the models as values.
    """
    # suffix = "" if suffix == "full" else suffix

    run_path = Path(run_path)

    model_paths = {
        model: list(run_path.glob(f"{prefix}{model}*{suffix}/{model_type}.pkl"))
        for model in model_names
    }

    return model_paths


# Metric Stuff
## Need this specifically correlation because native scipy pearsonr gives both the correlation and p-value which doesn't sit well when creating dataframes


def get_corr(preds: ndarray, ytest: ndarray, method: str) -> float:
    match method:
        case "pearson":
            return pearsonr(preds, ytest).statistic
        case "spearman":
            return spearmanr(preds, ytest).statistic
        case _:
            raise ValueError("Invalid method. Must be one of ['pearson', 'spearman']")
