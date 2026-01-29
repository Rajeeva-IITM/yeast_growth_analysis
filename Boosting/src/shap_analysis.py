from pathlib import Path
from typing import List

import hydra
import polars as pl
import shap
from dotenv import load_dotenv
from omegaconf import DictConfig
from rich import print
from sklearn.preprocessing import StandardScaler
from utils import get_model, get_model_paths

load_dotenv()


def get_shap_values(
    conf: DictConfig, model_path: Path, data_path: Path, fold: int
) -> pl.DataFrame:
    """Compute SHAP values for a given model and dataset.

    Parameters
    ----------
    model_path : Path
        Path to the pickle file containing the model.
    data_path : Path
        Path to the csv file containing the data.
    fold : int
        The number of the fold to compute SHAP values for.

    Returns
    -------
    pl.DataFrame
        A DataFrame containing the mean absolute SHAP values for each feature,
        sorted in descending order.
    """
    model = get_model(model_path)
    data: pl.DataFrame = pl.read_ipc(data_path)
    X = data.drop(["Condition", "Strain", "Phenotype"]).to_pandas()
    X_std = StandardScaler().fit_transform(X)
    y = data["Phenotype"].to_numpy()

    # instantiate explainer
    explainer: shap.TreeExplainer = hydra.utils.instantiate(
        conf.explainer, model, X_std
    )

    # compute shap values
    shap_values = explainer.shap_values(X_std, y, approximate=True)
    shap_df = pl.DataFrame(data=shap_values, schema=X.columns.to_list())
    shap_mean = (
        shap_df.with_columns(pl.all().abs().mean())
        .unique()
        .melt()
        .with_columns(Fold=pl.lit(fold))
        .rename({"value": "Value", "variable": "Feature"})
        .sort("Value", descending=True)
    )

    return shap_mean


def get_shap_folds(
    conf: DictConfig, model_paths: List[Path], data_path: Path
) -> pl.DataFrame:
    """Compute SHAP values for multiple models from the same dataset and concatenate them. Useful
    for when when you have models trained on different folds of the underlying data.

    Parameters
    ----------
    conf : DictConfig
        Configuration object containing the paths to the models and the data.
    model_paths : Path
        Paths to the pickle files containing the models.
    data_path : Path
        Path to the csv file containing the data.

    Returns
    -------
    pl.DataFrame
        A DataFrame containing the mean absolute SHAP values for each feature,
        sorted in descending order, computed for each model and concatenated.
    """
    df_folds = []
    for fold, model_path in enumerate(model_paths):
        print(f"Fold: {fold}")
        df_fold = get_shap_values(conf, model_path, data_path, fold)
        df_folds.append(df_fold)

    final_df: pl.DataFrame = pl.concat(df_folds)
    final_df = final_df.group_by("Feature", maintain_order=True).agg(
        pl.col("Value").mean()
    )
    return final_df


@hydra.main(  # pyrefly: ignore
    config_path="../configs/", version_base="1.3", config_name="interpret"
)
def main(conf: DictConfig) -> None:
    """Main entry point for the script.

    The script takes in a configuration, reads the models and data from the configuration,
    computes the SHAP values for each model, and saves the results to a CSV file.

    Parameters
    ----------
    conf : DictConfig
        The configuration object, passed in by Hydra.

    Returns
    -------
    None
    """

    print(conf)

    out_path = Path(conf.out_path)
    model_type = conf.model_load_keys.model_type
    model_paths = get_model_paths(**conf.model_load_keys)

    assert (
        len(conf.model_load_keys.model_names) <= len(conf.data_paths)
    ), f"Mismatch between data {len(conf.data_paths)} and model {len(conf.model_load_keys.model_names)} and model names"

    print(model_paths)

    for model_name, model_path in model_paths.items():
        shap_df = get_shap_folds(conf, model_path, conf.data_paths.get(model_name))
        shap_df.write_csv(out_path / f"shap_{model_name}_{model_type}.csv")


if __name__ == "__main__":
    main()
