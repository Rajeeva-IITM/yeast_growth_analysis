# Compare performance of any two model

from pathlib import Path
from typing import Union

import numpy as np
import hydra
import polars as pl

from dotenv import load_dotenv
from omegaconf import DictConfig
from rich import print

from utils import get_model, get_model_paths, get_data
from sklearn.preprocessing import StandardScaler

load_dotenv()


def get_preds(
    path_of_model: Union[str, Path],
    data_path: Union[str, Path],
    run_type: str,
    regression: bool,
) -> pl.DataFrame:
    """A function that generates predictions using a trained model.

    Parameters:
        model_path (str): Path to the trained model file.
        data_path (str): Path to the data file for prediction.
        run_type (str): One of `full`, `geno_only`, `chem_only`
        regression (bool): Whether the model is regression or classification

    Returns:
        pl.DataFrame: A DataFrame containing the actual Phenotype, predicted Phenotype, and Condition.
    """
    # print("We got here!")
    model = get_model(path_of_model)
    data = pl.read_ipc(data_path)

    match run_type:
        case "full" | "geno_only" | "chem_only":
            X, y = get_data(data_path=data_path, run_type=run_type, return_as_Xy=True)
        case "dummy":
            X, y = get_data(
                data_path=data_path, run_type="full", return_as_Xy=True
            )  # dummy is the same as full, because it doesn't care about features

    X = StandardScaler().fit_transform(X)
    preds = model.predict(X)

    if not regression:
        preds = np.where(preds > 0.5, 1, 0)

    return pl.DataFrame(
        {"Phenotype": y, "Preds": preds, "Condition": data["Condition"]}
    )


def get_preds_kfold(
    paths_of_model: Union[str, Path.glob], data_path, run_type, regression: bool
) -> pl.DataFrame:
    """A function that generates predictions using a trained model.

    Parameters:
        paths_of_model (str): Path to the trained model file.
        data_path (str): Path to the data file for prediction.
        run_type (str): One of `full`, `geno_only`, `chem_only`
        regression (bool): Whether the model is regression or classification

    Returns:
        pl.DataFrame: A DataFrame containing the actual Phenotype, predicted Phenotype, and Condition.
    """
    preds = []
    for i, path_of_model in enumerate(paths_of_model):
        result = get_preds(path_of_model, data_path, run_type, regression)
        result = result.with_columns(pl.lit(i).alias("Fold"))
        preds.append(result)
    return pl.concat(preds)


def eval_model(conf: DictConfig, pred_df: pl.DataFrame) -> pl.DataFrame:
    """Evaluates the performance of a model based on its predictions.

    Parameters:
        conf (DictConfig): A configuration object with the relevant metrics
        pred_df (pl.DataFrame): A DataFrame containing the predicted values,
            actual values, and conditions.

    Returns:
        pl.DataFrame: A DataFrame containing the accuracy, f1 score, auc roc score,
            and mathews correlation coefficient for each condition and fold.
    """

    result_df = pred_df.group_by("Condition", "Fold", maintain_order=True).map_groups(
        lambda x: pl.DataFrame(
            {
                "Condition": x["Condition"].unique(maintain_order=True),
                "Fold": x["Fold"].unique(maintain_order=True),
                # "accuracy": accuracy_score(x["Phenotype"], x["Preds"]),
                # "f1": f1_score(x["Phenotype"], x["Preds"]),
                # "auc_roc": roc_auc_score(x["Phenotype"], x["Preds"]),
                # "mathews": matthews_corrcoef(x["Phenotype"], x["Preds"]),
            }
            | {
                metric: hydra.utils.call(
                    conf.metrics.get(metric), x["Phenotype"], x["Preds"]
                )
                for metric in conf.metrics
            }
        )
    )

    return result_df


def get_results(conf: DictConfig):
    """Runs the models in the `run_path` directory on the data in the `data_paths` dictionary, and
    saves the results to the `out_path` directory.

    Parameters:
        conf (DictConfig): A configuration object with the following required
            keys:
                - run_path (Path): The path to the directory containing the models.
                - data_paths (Dict[str, Path]): A dictionary mapping data names to
                    paths to the data files.
                - out_path (Path): The path to the directory to save the results to.
                - model_types (List[str]): A list of model types to run.
                - run_type (str): The suffix of the model name to distinguish between
                    different runs (e.g. 'full', 'geno_only', 'chem_only', 'dummy').
                - regression (bool): Whether the model is regression or classification.

    Returns:
        None
    """
    # run_path = Path(conf.run_path)
    # data_path = Path(conf.data_paths)
    out_path = Path(conf.out_path)
    model_type = conf.model_load_keys.model_type
    model_paths = get_model_paths(**conf.model_load_keys)

    assert (
        len(conf.model_load_keys.model_names) <= len(conf.data_paths)
    ), f"Mismatch between data {len(conf.data_paths)} and model {len(conf.model_load_keys.model_names)} and model names"

    print(model_paths)

    for model_name, model_path in model_paths.items():
        for data_name, data_path in conf.data_paths.items():
            print(model_name, data_name)
            result_df = get_preds_kfold(
                model_path, data_path, conf.run_type, conf.regression
            )
            metric_df = eval_model(conf, result_df)

            result_df.write_parquet(
                out_path / f"predictions_{model_name}_{data_name}_{model_type}.parquet"
            )

            metric_df.write_csv(
                out_path / f"metrics_{model_name}_{data_name}_{model_type}.csv"
            )


@hydra.main(config_path="../configs/", version_base="1.3", config_name="eval")
def main(conf: DictConfig) -> None:
    """The main entry point of the script.

    Args:
        conf (DictConfig): The configuration object, passed in by Hydra.

    Returns:
        None
    """
    print(conf)
    get_results(conf)


if __name__ == "__main__":
    main()
