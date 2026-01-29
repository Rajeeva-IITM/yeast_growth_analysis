# Want to check what kind of data when sampled creates optimal performance

import hydra
import numpy as np
import polars as pl
import polars.selectors as cs
import optuna
import pickle

from warnings import warn

from pathlib import Path
from typing import Tuple, Union
from dotenv import load_dotenv
from omegaconf import DictConfig
import lightgbm as lgb

# from sklearn.metrics import roc_auc_score
# from sklearn.model_selection import KFold, train_test_split
from sklearn.preprocessing import StandardScaler
from rich.console import Console

from scipy.spatial.distance import pdist, squareform

load_dotenv()

console = Console(record=True)


def data_sampler(
    # config: DictConfig,
    df: pl.DataFrame,
    n_samples: int,
    grouping: bool = False,
    grouping_strategy: str = None,
    num_groups: Union[float, int] = 0.1,
    seed: int = 42,
) -> Tuple[pl.DataFrame, np.array]:
    """Samples data from a DataFrame based on the specified criteria.

    Args:
        df (pl.DataFrame): The DataFrame to sample from.
        n_samples (int): The number of samples to extract.
        grouping (bool, optional): Whether to group the samples. Defaults to False.
        grouping_strategy (str, optional): The strategy for grouping the samples. Defaults to None.
        num_groups (Union[float, int], optional): The number of groups to consider. Defaults to 0.1.
        seed (int, optional): The seed for random number generation. Defaults to 42.

    Returns:
        Tuple[pl.DataFrame, np.array]: A tuple containing the sampled data and the corresponding labels.

    Raises:
        AssertionError: If the number of samples is invalid.
        AssertionError: If the number of groups is invalid.
        ValueError: If the grouping_strategy is invalid.
    """
    assert (
        n_samples > 0 and isinstance(n_samples, int) and n_samples <= len(df)
    ), "Invalid number of samples"

    assert isinstance(num_groups, (float, int)), "`num_groups` must be float or int"

    # if grouping_strategy is not None:
    #     assert (
    #         grouping
    #     ), "Chemical aware strategy can only be used when grouping is True"

    if not grouping:
        samples = df.sample(n=n_samples, seed=seed)
        X, y = samples.drop(["Condition", "Strain", "Phenotype"]), samples["Phenotype"]

    else:
        match grouping_strategy:
            case "condition":
                unique_conditions = df["Condition"].value_counts()

                if isinstance(num_groups, int):
                    assert (
                        num_groups > 0 and num_groups < len(unique_conditions)
                    ), "Number of groups must be between 0 and the number of unique conditions"

                else:
                    assert (
                        0 < num_groups < 1
                    ), "Number of groups must be between 0 and 1"
                    num_groups = np.ceil(len(unique_conditions) * num_groups)
                    # This is not ideal because sometimes the requested number of groups may not be enough for the number of samples
                    # In this case, I'm sampling with replacement

                conditions_considered = unique_conditions.sample(num_groups, seed=seed)[
                    "Condition"
                ].to_list()

                if (
                    unique_conditions.filter(
                        pl.col("Condition").is_in(conditions_considered)
                    ).sum()["count"][0]
                    < n_samples
                ):
                    warn(
                        f"Sampling with replacement is done. The number of groups {num_groups} is not enough for the number of samples {n_samples}"
                    )
                    samples = df.filter(
                        pl.col("Condition").is_in(conditions_considered)
                    ).sample(n=n_samples, seed=seed, with_replacement=True)

                else:
                    samples = df.filter(
                        pl.col("Condition").is_in(conditions_considered)
                    ).sample(n=n_samples, seed=seed)

                X, y = (
                    samples.drop(["Condition", "Strain", "Phenotype"]),
                    samples["Phenotype"],
                )

            case "strain":
                unique_strains = df["Strain"].value_counts()

                if isinstance(num_groups, int):
                    assert (
                        num_groups > 0 and num_groups < len(unique_strains)
                    ), "Number of groups must be between 0 and the number of unique strains"

                else:
                    assert (
                        0 < num_groups < 1
                    ), "Number of groups must be between 0 and 1"
                    num_groups = np.ceil(len(unique_strains) * num_groups)
                    # This is not ideal because sometimes the requested number of groups may not be enough for the number of samples
                    # In this case, I'm sampling with replacement

                strains_considered = unique_strains.sample(num_groups, seed=seed)[
                    "Strain"
                ].to_list()

                if (
                    unique_strains.filter(
                        pl.col("Strain").is_in(strains_considered)
                    ).sum()["count"][0]
                    < n_samples
                ):
                    warn(
                        f"Sampling with replacement is done. The number of groups {num_groups} is not enough for the number of samples {n_samples}"
                    )
                    samples = df.filter(
                        pl.col("Strain").is_in(strains_considered)
                    ).sample(n=n_samples, seed=seed, with_replacement=True)

                else:
                    samples = df.filter(
                        pl.col("Strain").is_in(strains_considered)
                    ).sample(n=n_samples, seed=seed)

                X, y = (
                    samples.drop(["Condition", "Strain", "Phenotype"]),
                    samples["Phenotype"],
                )

            case "intelligent_strain":
                # A little more involved
                # Obtain samples from strains that are most representative
                # print("Sampling from representative strains")
                unique_strains = df["Strain"].value_counts()
                # strains = df["Strain"].unique(maintain_order=True).to_numpy()
                strain_data = (
                    df.select(cs.starts_with("Y"))
                    # .unique(maintain_order=True)
                    .to_numpy()
                )
                avg_dist_btw_strains = squareform(
                    pdist(strain_data, metric="cityblock")
                ).mean(axis=0)

                representative_strains = np.argsort(avg_dist_btw_strains)[
                    : int(
                        num_groups * len(avg_dist_btw_strains)
                        if isinstance(num_groups, float)
                        else num_groups
                    )
                ]

                if (
                    unique_strains.filter(
                        pl.col("Strain").is_in(representative_strains)
                    ).sum()["count"][0]
                    < n_samples
                ):
                    warn(
                        f"Sampling with replacement is done. The number of groups {num_groups} is not enough for the number of samples {n_samples}"
                    )
                    samples = df.filter(
                        pl.col("Strain").is_in(representative_strains)
                    ).sample(n=n_samples, seed=seed, with_replacement=True)
                else:
                    samples = df.filter(
                        pl.col("Strain").is_in(representative_strains)
                    ).sample(n=n_samples, seed=seed)

                X, y = (
                    samples.drop(["Condition", "Strain", "Phenotype"]),
                    samples["Phenotype"],
                )
            case _:
                ValueError("Invalid grouping_strategy")

    return X, y.to_numpy()


def tune_sampler(
    trial: optuna.Trial,
    train_df: pl.DataFrame,
    test_df: pl.DataFrame,
    config: DictConfig,
):
    """Tunes the data_sampler parameters using Optuna's trial object.

    Parameters
    ----------
    trial : optuna.Trial
        The Optuna trial object.
    train_df : pl.DataFrame
        The input DataFrame.
    test_df : pl.DataFrame
        The test DataFrame.
    config : DictConfig
        The configuration dictionary.

    Returns
    -------
    n_samples : int
        The optimal number of samples.
    mean_score : float
        The mean score over the cross-validation folds.
    """
    X_test, y_test = (
        test_df.drop(["Condition", "Strain", "Phenotype"]).to_numpy(),
        test_df["Phenotype"].to_numpy(),
    )
    X_test = StandardScaler().fit_transform(X_test)
    val_dataset = lgb.Dataset(X_test, y_test, free_raw_data=False)

    max_samples = train_df.shape[0]  # Max should be # of the training samples

    sampler_params = {
        "df": train_df,
        "n_samples": trial.suggest_int(
            "n_samples",
            10,
            max_samples,
        ),
        "grouping": trial.suggest_categorical("grouping", [True, False]),
        "grouping_strategy": trial.suggest_categorical(
            "grouping_strategy",
            [
                "condition",
                "strain",
            ],
        ),
        "num_groups": trial.suggest_float("num_groups", 0.1, 0.9),
        # "seed": config.seed,
    }

    scores = []

    for i in range(config.repeat_num):
        sampler_params["seed"] = i if config.repeat_samples else config.seed

        sampled_X, sampled_y = data_sampler(**sampler_params)
        sampled_X = StandardScaler().fit_transform(sampled_X.to_numpy())

        train_dataset = lgb.Dataset(sampled_X, sampled_y, free_raw_data=False)

        model_params = {
            "objective": config.model_params.objective,
            "verbose": -10,
            "early_stopping_rounds": 10,
            "lambda_l1": config.model_params.lambda_l1,
            "lambda_l2": config.model_params.lambda_l2,
            "num_leaves": config.model_params.num_leaves,
            "feature_fraction": config.model_params.feature_fraction,
            "bagging_fraction": config.model_params.bagging_fraction,
            "bagging_freq": config.model_params.bagging_freq,
            "min_child_samples": config.model_params.min_child_samples,
            "learning_rate": config.model_params.learning_rate,
            "metrics": config.model_params.metrics,
            "random_seed": i
            if not config.model_params.seed
            else config.model_params.seed,  # Change seed every iteration
        }

        n_estimators = config.model_params.n_estimators

        model = lgb.train(
            params=model_params,
            train_set=train_dataset,
            num_boost_round=n_estimators,
            valid_sets=[val_dataset],
            # valid_names=["val"],
            # feval=eval_metric
        )

        # console.print(model.best_score["valid_0"])

        score = model.best_score["valid_0"][config.model_params.metrics]

        scores.append(score)

    return sampler_params["n_samples"], np.mean(scores)


def verify_path(path: str):
    """Checks if a path exists and creates it if it doesn't.

    Args:
        path (str): The path to check/create.

    Returns:
        None
    """
    path = Path(path)
    if not path.exists():
        path.mkdir(parents=True)
        console.print(f"Created directory at [red]{path}[/red]", justify="center")


@hydra.main(version_base=None, config_path="../configs/", config_name="sampler_conf")
def main(conf: DictConfig):
    """Main function to run the hyperparameter tuning for the sampler.

    Parameters:
        conf (DictConfig): The configuration from the YAML file.

    Returns:
        None
    """
    console.log("Processing data", style="bold red", justify="center")

    verify_path(conf.data.savedir)  # Create directory if it doesn't exist

    train_df = pl.read_ipc(conf.data.train_data)
    test_df = pl.read_ipc(conf.data.test_data)

    storage_path = conf.data.savedir + "tune_study_journal.txt"
    storage = optuna.storages.JournalStorage(
        optuna.storages.JournalFileStorage(storage_path)
    )

    sampler_study = optuna.create_study(
        storage=storage,
        study_name="sampler_study",
        directions=conf.study_directions,
        pruner=optuna.pruners.HyperbandPruner(),
        sampler=optuna.samplers.TPESampler(seed=conf.seed),
        load_if_exists=True,
    )

    sampler_study.optimize(
        lambda trial: tune_sampler(trial, train_df, test_df, conf),
        n_trials=conf.n_trials,
        # catch=(ValueError)
    )

    with open(conf.data.savedir + "sampler_study.pkl", "wb") as f:
        pickle.dump(sampler_study.best_trials, f)
    # best_params = sampler_study.best_trials[0].params
    # print(best_params)


if __name__ == "__main__":
    # # Write a test for the data_sampler

    # df = pl.read_ipc(
    #     "/home/rajeeva/Project/data/chem_subsets/new_latent/carbons/carbons_bloom2013_clf_3_pubchem.feather"
    # )
    # n_samples = 100
    # X, y = data_sampler(df, n_samples, grouping=True, grouping_strategy="strain")
    # print(X)
    # print(y)

    main()
