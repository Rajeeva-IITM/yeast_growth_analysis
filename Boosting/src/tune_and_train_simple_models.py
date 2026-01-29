# Running the chemical subsets results again because consistency and I'm an idiot
import logging
import pickle
from pathlib import Path
from typing import Callable

import hydra
import numpy as np
import optuna
import polars as pl
import polars.selectors as cs

from sklearnex import patch_sklearn

patch_sklearn()

from sklearn.neighbors import KNeighborsRegressor  # noqa: E402
from sklearn.ensemble import RandomForestClassifier  # noqa: E402

# from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, ElasticNet  # noqa: E402
from sklearn.svm import SVC, SVR  # noqa: E402
from dotenv import load_dotenv  # noqa: E402
from omegaconf import DictConfig  # noqa: E402
from rich.console import Console  # noqa: E402
from rich.logging import RichHandler  # noqa: E402

# from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold, train_test_split  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

load_dotenv()

console = Console(record=True)

logger = logging.getLogger()
logger.addHandler(RichHandler(console=console))

optuna.logging.disable_default_handler()  # Stop showing logs in sys.stderr.
optuna.logging.enable_propagation()  # Propagate logs to the root logger.

intel_logger = logging.getLogger("sklearnex")
# intel_logger.removeHandler(intel_logger.handlers[0])
# all_logs = logging.FileHandler('/data/rajeeva/Boosting-yeast_growth_pred/outputs/sklearnex.log')
# all_logs.setFormatter(logging.Formatter('%(message)s'))
intel_logger.setLevel(logging.WARNING)
# intel_logger.addHandler(all_logs)


def tune_LogRegression(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the Logistic Regression model using Optuna's trial object.

    Args:
        trial (optuna.Trial): The Optuna trial object.
        X (numpy.ndarray | DataFrame): The input features.
        y (numpy.ndarray | DataFrame): The target variable.
        config (DictConfig): The configuration dictionary.

    Returns:
        float: The average score over the cross-validation folds.
    """
    penalty = trial.suggest_categorical(
        "penalty",
        [
            "l1",
            "l2",
            "elasticnet",
        ],
    )
    C = trial.suggest_float("C", 1e-3, 1e3)

    params = {
        "penalty": penalty,
        "C": C,
        "max_iter": 10000,
    }

    if penalty == "elasticnet":
        l1_ratio = trial.suggest_float("l1_ratio", 0.0, 1.0)
        params["l1_ratio"] = l1_ratio

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )
    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = LogisticRegression(**params)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        scores.append(hydra.utils.call(config.metric, _args_=(y_test, y_pred)))

    return np.mean(scores)


def tune_SVM(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the Support Vector Machine (SVM) model using Optuna's trial object.

    Args:
        trial (optuna.Trial): The Optuna trial object.
        X (numpy.ndarray | DataFrame): The input features.
        y (numpy.ndarray | DataFrame): The target variable.
        config (DictConfig): The configuration dictionary.

    Returns:
        float: The average score over the cross-validation folds.
    """
    params = {
        "C": trial.suggest_float("C", 1e-3, 1e3),
        "kernel": trial.suggest_categorical(
            "kernel", ["linear", "poly", "rbf", "sigmoid"]
        ),
        "gamma": trial.suggest_categorical("gamma", ["scale", "auto"]),
    }

    if params["kernel"] == "poly":
        degree = trial.suggest_int("degree", 2, 7)
        params["degree"] = degree

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )
    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = SVC(**params)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        scores.append(hydra.utils.call(config.metric, _args_=(y_test, y_pred)))

    return np.mean(scores)


def tune_RF(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the Random Forest Classifier model using Optuna's trial object.

    Args:
        trial (optuna.Trial): The Optuna trial object.
        X (numpy.ndarray | DataFrame): The input features.
        y (numpy.ndarray | DataFrame): The target variable.
        config (DictConfig): The configuration dictionary.

    Returns:
        float: The average score over the cross-validation folds.
    """
    n_estimators = trial.suggest_int("n_estimators", 10, 1000, log=True)
    criterion = trial.suggest_categorical("split_criterion", ["gini", "entropy"])
    max_samples = trial.suggest_float("max_samples", 0.1, 1.0)
    max_features = trial.suggest_float("max_features", 0.1, 1.0)
    max_depth = trial.suggest_int("max_depth", 2, 12)
    min_samples_leaf = trial.suggest_int("min_samples_leaf", 5, 100, log=True)
    # ccp_alpha = trial.suggest_float("ccp_alpha", 1e-7, 1.0 , log=True)

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )
    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = RandomForestClassifier(
            n_estimators=n_estimators,
            criterion=criterion,
            max_samples=max_samples,
            max_features=max_features,
            max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            # ccp_alpha=ccp_alpha,
            # n_streams=16,
            random_state=config.seed,
        )
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        scores.append(hydra.utils.call(config.metric, _args_=(y_test, y_pred)))

    return np.mean(scores)


# Regression functions


def tune_ElasticNet(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the Elastic Net regression model using Optuna's trial object.

    Args:
        trial (optuna.Trial): The Optuna trial object.
        X (numpy.ndarray | DataFrame): The input features.
        y (numpy.ndarray | DataFrame): The target variable.
        config (DictConfig): The configuration dictionary.

    Returns:
        float: The average score over the cross-validation folds.
    """
    params = {
        "alpha": trial.suggest_float("alpha", 1e-3, 10, log=True),
        "l1_ratio": trial.suggest_float("l1_ratio", 0.0, 1.0),
        # "n_jobs":
    }

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )
    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = ElasticNet(**params)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        scores.append(hydra.utils.call(config.metric, _args_=(y_test, y_pred)))

    return np.mean(scores)


def tune_SVR(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the Support Vector Regression (SVR) model using Optuna's trial object.

    Args:
        trial (optuna.Trial): The Optuna trial object.
        X (numpy.ndarray | DataFrame): The input features.
        y (numpy.ndarray | DataFrame): The target variable.
        config (DictConfig): The configuration dictionary.

    Returns:
        float: The average score over the cross-validation folds.
    """
    params = {
        "C": trial.suggest_float("C", 1e-3, 1e3),
        "kernel": trial.suggest_categorical(
            "kernel",
            [
                "poly",
                "rbf",
            ],
        ),
        "gamma": trial.suggest_categorical("gamma", ["scale", "auto"]),
    }

    if params["kernel"] == "poly":
        degree = trial.suggest_int("degree", 2, 7)
        params["degree"] = degree

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )
    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = SVR(**params)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        scores.append(hydra.utils.call(config.metric, _args_=(y_test, y_pred)))

    return np.mean(scores)


def tune_Neighbours(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the K-Neighbors regression model using Optuna's trial object.

    Args:
        trial (optuna.Trial): The Optuna trial object.
        X (numpy.ndarray | DataFrame): The input features.
        y (numpy.ndarray | DataFrame): The target variable.
        config (DictConfig): The configuration dictionary.

    Returns:
        float: Mean score over the cross-validation folds.
    """
    params = {
        "n_neighbors": trial.suggest_int("n_neighbors", 5, 500),
        "weights": trial.suggest_categorical("weights", ["uniform", "distance"]),
        "metric": "euclidean",
    }

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )

    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = KNeighborsRegressor(**params)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        scores.append(hydra.utils.call(config.metric, _args_=(y_test, y_pred)))


def run_study(
    study_name: str, tune_objective: Callable, conf: DictConfig, X, y, n_jobs: int = 1
):
    """Runs an Optuna study to optimize hyperparameters for a given model.

    Args:
        study_name (str): The name of the study.
        tune_objective (Callable): The function to optimize. Should take in a
            trial object, the training data, and the configuration object.
        conf (DictConfig): The configuration for the study.
        X: The training features.
        y: The training target variable.
        n_jobs (int, optional): The number of jobs to run in parallel. Defaults to 1.

    Returns:
        optuna.study.Study: The study object.
    """
    pruner = optuna.pruners.HyperbandPruner()
    sampler = optuna.samplers.TPESampler(seed=conf.seed)
    study = optuna.create_study(
        study_name=study_name,
        direction="maximize",
        pruner=pruner,
        sampler=sampler,
    )
    study.optimize(
        lambda trial: tune_objective(trial, X, y, conf),
        n_trials=conf.n_trials,
        catch=(ValueError),
        n_jobs=n_jobs,  # Carefully use this parameter. Avoid when using multithreaded algorithms
    )

    return study


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


@hydra.main(version_base=None, config_path="../configs/", config_name="conf")
def main(conf: DictConfig):
    """The main function that serves as the entry point for the program.

    Args:
        conf (DictConfig): The configuration object containing various settings.

    Returns:
        None
    """

    # Process Data

    console.log("Processing data", style="bold red", justify="center")

    verify_path(conf.data.savedir)

    df = pl.read_ipc(conf.data.path)

    match conf.run_type:
        case "geno_only":
            df = df.select(
                pl.col("Strain"),
                cs.starts_with("Y"),  # Only the genotype columns
                pl.col("Condition"),
                pl.col("Phenotype"),
            )
        case "chem_only":
            df = df.select(
                pl.col("Strain"),
                cs.contains("latent"),  # Only the chemical information
                pl.col("Condition"),
                pl.col("Phenotype"),
            )
        case "full":
            pass
        case _:
            raise ValueError(
                "Invalid run_type. Must be one of ['geno_only', 'chem_only', 'full']"
            )

    X = (
        df.drop(
            ["Phenotype", "Condition", "Strain"],
        )
        .with_columns(pl.all().cast(pl.Float32))
        .to_numpy()
    )
    y = df["Phenotype"].to_numpy()

    Xtrain, Xtest, ytrain, ytest = train_test_split(
        X, y, test_size=0.2, random_state=conf.kfold_params.seed
    )

    scaler = StandardScaler()
    Xtrain = scaler.fit_transform(Xtrain)
    Xtest = scaler.transform(Xtest)

    console.log("Training models", style="bold green", justify="center")

    if conf.regression:
        # Optimizing elastic Net

        elasticNet_study = run_study(
            study_name=conf.models.ElasticNet.name,
            tune_objective=tune_ElasticNet,
            conf=conf,
            X=Xtrain,
            y=ytrain,
            n_jobs=5,  # Carefully use this parameter. Avoid when using multithreaded algorithms
        )

        # Train best model
        best_model = ElasticNet(**elasticNet_study.best_trial.params)
        best_model.fit(Xtrain, ytrain)

        # verify_path(conf.models.ElasticNet.model_savename)
        pickle.dump(best_model, open(conf.models.ElasticNet.model_savename, "wb"))

        console.log("Elastic Net Done", style="bold green", justify="center")

        # # Optimizing N Neighbours

        # n_neighbours_study = run_study(
        #     study_name=conf.models.NeighboursR.model_savename,
        #     tune_objective=tune_Neighbours,
        #     conf=conf,
        #     X=Xtrain,
        #     y=ytrain,
        #     n_jobs = 25 # Carefully use this parameter. Avoid when using multithreaded algorithms
        # )

        # best_model = KNeighborsRegressor(**n_neighbours_study.best_trial.params)
        # best_model.fit(Xtrain, ytrain)

        # pickle.dump(best_model, open(conf.models.NeighboursR.model_savename, "wb"))

        # console.log("K-Neighbors Done", style="bold green", justify="center")

        # Optimizing SVR - Takes too long due to the large amount of data samples

        svr_study = run_study(
            study_name=conf.models.SVR.name,
            tune_objective=tune_SVR,
            conf=conf,
            X=Xtrain,
            y=ytrain,
            n_jobs=5,
        )

        # Train best model
        best_model = SVR(**svr_study.best_trial.params)
        best_model.fit(Xtrain, ytrain)

        # verify_path(conf.models.SVR.model_savename)
        pickle.dump(best_model, open(conf.models.SVR.model_savename, "wb"))

        console.log("SVR Done", style="bold green", justify="center")

    else:
        # Optimizing for Logistic Regression
        logReg_study = run_study(
            study_name=conf.models.LogReg.name,
            tune_objective=tune_LogRegression,
            conf=conf,
            X=Xtrain,
            y=ytrain,
        )

        # Train best model
        best_model = LogisticRegression(**logReg_study.best_trial.params)
        best_model.fit(Xtrain, ytrain)

        # verify_path(conf.models.LogReg.model_savename)
        pickle.dump(best_model, open(conf.models.LogReg.model_savename, "wb"))

        console.log("Logistic Regression Done", style="bold green", justify="center")

        # SVM
        svm_study = run_study(
            study_name=conf.models.SVM.name,
            tune_objective=tune_SVM,
            conf=conf,
            X=Xtrain,
            y=ytrain,
        )
        best_model = SVC(**svm_study.best_trial.params)
        best_model.fit(Xtrain, ytrain)

    # verify_path(conf.models.SVM.model_savename)
    pickle.dump(best_model, open(conf.models.SVM.model_savename, "wb"))

    console.log("SVM Done", style="bold green", justify="center")

    # Random Forest - This was causing a lot of issues
    # rf_study = run_study(
    #     study_name=conf.models.RandomForest.name,
    #     tune_objective=tune_RF,
    #     conf=conf,
    #     X=Xtrain,
    #     y=ytrain,
    # )
    # best_model = RandomForestClassifier(**rf_study.best_trial.params)
    # best_model.fit(Xtrain, ytrain)

    # # verify_path(conf.models.RandomForest.model_savename)
    # pickle.dump(best_model, open(conf.models.RandomForest.model_savename, "wb"))
    # console.log("Random Forest Done", style="bold green", justify="center")

    console.save_html(conf.data.savedir + "/run_report.html")


if __name__ == "__main__":
    main()
