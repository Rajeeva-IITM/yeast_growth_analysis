# Running the chemical subsets results again because consistency and I'm an idiot
import logging
import pickle
import sys
from pathlib import Path
from typing import Callable

import hydra
import lightgbm as lgb
import numpy as np
import optuna
from dotenv import load_dotenv
from omegaconf import DictConfig
from rich.console import Console
from rich.logging import RichHandler
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import KFold, train_test_split
from sklearn.preprocessing import StandardScaler
from train import train_booster
from utils import get_data

load_dotenv()

console = Console(record=True)

logger = logging.getLogger()
logger.addHandler(RichHandler(console=console))

optuna.logging.disable_default_handler()  # Stop showing logs in sys.stderr.
optuna.logging.enable_propagation()  # Propagate logs to the root logger.


def create_dummy(
    strategy: str, Xtrain, ytrain, config: DictConfig, Xtest=None, ytest=None
):
    # splits = KFold(
    #     n_splits=config.kfold_params.n_splits,
    #     shuffle=config.kfold_params.shuffle,
    #     random_state=config.kfold_params.seed,
    # )

    """Creates a dummy classifier given the strategy and fits it to the training data. If test data
    is provided, it will also evaluate the model on the test data and log the score.

    Parameters
    ----------
    strategy : str
        The strategy for the dummy classifier. Can be one of 'most_frequent',
        'uniform', 'constant', or 'stratified'.
    Xtrain : array-like of shape (n_samples, n_features)
        The training data.
    ytrain : array-like of shape (n_samples,)
        The target values for the training data.
    config : DictConfig
        The configuration dictionary.
    Xtest : array-like of shape (n_samples, n_features), optional
        The test data.
    ytest : array-like of shape (n_samples,), optional
        The target values for the test data.

    Returns
    -------
    model : DummyClassifier
        The trained dummy classifier.
    """
    if (Xtest is not None) and (ytest is not None):
        # test_scores = []
        X_test, y_test = Xtest, ytest

    model = DummyClassifier(strategy=strategy, random_state=config.seed)
    model.fit(Xtrain, ytrain)

    y_pred_test = model.predict(X_test)
    score = hydra.utils.call(config.metric, _args_=(y_test, y_pred_test))

    console.log(
        f"Dummy classifier {strategy} score: {score:.2f}",
        style="bold red",
        justify="center",
    )

    return model


# A function to tune a LightGBM model
def tune_LGBM(
    trial: optuna.Trial, Xtrain, ytrain, config: DictConfig, Xtest=None, ytest=None
):
    """Tunes the LightGBM model using Optuna's trial object.

    Args:
        trial (optuna.Trial): The Optuna trial object.
        Xtrain (numpy.ndarray | DataFrame): The input features.
        ytrain (numpy.ndarray | DataFrame): The target variable.
        config (DictConfig): The configuration dictionary.
        Xtest (numpy.ndarray | DataFrame, optional): The test input features. Defaults to None.
        ytest (numpy.ndarray | DataFrame, optional): The test target variable. Defaults to None.

    Returns:
        float: The average ROC AUC score over the cross-validation folds.
    """
    # Boosting parameters
    params = {
        "objective": config.model_params.objective,
        "verbosity": config.model_params.verbosity,
        "force_col_wise": config.model_params.force_col_wise,
        "early_stopping_rounds": config.model_params.early_stopping_rounds,
        "num_threads": config.model_params.num_threads,
        "boosting_type": config.model_params.boosting_type,
        "device_type": config.model_params.device_type,
        "gpu_use_dp": config.model_params.gpu_use_dp,
        "seed": config.seed,
        "lambda_l1": trial.suggest_float("lambda_l1", **config.model_params.lambda_l1),
        "lambda_l2": trial.suggest_float("lambda_l2", **config.model_params.lambda_l2),
        "num_leaves": trial.suggest_int("num_leaves", **config.model_params.num_leaves),
        "feature_fraction": trial.suggest_float(
            "feature_fraction", **config.model_params.feature_fraction
        ),
        "bagging_fraction": trial.suggest_float(
            "bagging_fraction", **config.model_params.bagging_fraction
        ),
        "bagging_freq": trial.suggest_int(
            "bagging_freq", **config.model_params.bagging_freq
        ),
        "min_child_samples": trial.suggest_int(
            "min_child_samples", **config.model_params.min_child_samples
        ),
        "learning_rate": trial.suggest_float(
            "learning_rate", **config.model_params.learning_rate
        ),
    }

    n_estimators = trial.suggest_int("n_estimators", **config.model_params.n_estimators)

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )

    scores = []  # List of scores
    if (Xtest is not None) and (ytest is not None):
        # test_scores = []
        X_test, y_test = Xtest, ytest
        # test_dataset = lgb.Dataset(X_test, label=y_test)

    for train_idx, val_idx in splits.split(Xtrain, ytrain):
        X_train, y_train = Xtrain[train_idx], ytrain[train_idx]
        # console.print(Xtrain.dtype, X_test.dtype)
        # if (Xtest is  None) and (ytest is None):
        X_val, y_val = Xtrain[val_idx], ytrain[val_idx]
        val_dataset = lgb.Dataset(X_val, label=y_val)
        train_dataset = lgb.Dataset(X_train, label=y_train)
        model = lgb.train(
            params,
            train_dataset,
            num_boost_round=n_estimators,
            valid_sets=[val_dataset],
            # verbose_eval=False
        )

        y_pred_val = model.predict(
            X_val,
        )
        scores.append(
            hydra.utils.call(config.metric, _args_=(y_val, y_pred_val))
        )  # Include score for validation set

        # if separate test set is given
        if (Xtest is not None) and (ytest is not None):
            y_pred_test = model.predict(X_test)
            scores.append(
                hydra.utils.call(config.metric, _args_=(y_test, y_pred_test))
            )  # Include score for test set

    return np.mean(scores)


def run_study(
    study_name: str,
    tune_objective: Callable,
    conf: DictConfig,
    Xtrain,
    ytrain,
    Xtest=None,
    ytest=None,
):
    """Runs an Optuna study to optimize hyperparameters for a given model.

    Args:
        study_name (str): The name of the study.
        tune_objective (Callable): The function to optimize. Should take in a
            trial object, the training data, and the configuration object.
        conf (DictConfig): The configuration for the study.
        Xtrain: The training features.
        ytrain: The training target variable.
        Xtest (optional): The test features. Defaults to None.
        ytest (optional): The test target variable. Defaults to None.

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
        lambda trial: tune_objective(trial, Xtrain, ytrain, conf, Xtest, ytest),
        n_trials=conf.n_trials,
        # catch=(ValueError)
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

    Xtrain, ytrain = get_data(conf.data.path, run_type=conf.run_type, return_as_Xy=True)

    # If no separate test set is given
    if conf.testing.test_dataset is None:
        Xtrain, Xtest, ytrain, ytest = train_test_split(
            Xtrain,
            ytrain,
            test_size=conf.testing.test_frac,
            random_state=conf.seed,
        )

        scaler = StandardScaler()
        Xtrain = scaler.fit_transform(Xtrain)
        Xtest = scaler.transform(Xtest)

    else:
        Xtest, ytest = get_data(
            conf.testing.test_dataset, run_type=conf.run_type, return_as_Xy=True
        )
        scaler = StandardScaler()
        Xtrain = scaler.fit_transform(Xtrain)
        Xtest = scaler.transform(Xtest)

    console.log("Data processed", style="bold green", justify="center")

    if conf.dummy.run:
        console.log("Creating DummyClassifier", style="bold red", justify="center")

        model = create_dummy(
            strategy=conf.dummy.strategy,
            Xtrain=Xtrain,
            ytrain=ytrain,
            config=conf,
            Xtest=Xtest,
            ytest=ytest,
        )

        savename = "/" + (conf.dummy.strategy) + ".pkl"

        with open(conf.data.savedir + savename, "wb") as f:
            pickle.dump(model, f)

        console.save_html(conf.data.savedir + f"/{savename}_run_report.html")

        return None  # End the program
    else:
        pass

    console.log("Training models", style="bold green", justify="center")

    study_name = (
        conf.models.Boosting.name
        if conf.model_params.boosting_type == "gbdt"
        else conf.models.RandomForest.name
    )

    # Boosting
    boost_study = run_study(
        study_name=study_name,
        tune_objective=tune_LGBM,
        conf=conf,
        Xtrain=Xtrain,
        ytrain=ytrain,
        Xtest=Xtest,
        ytest=ytest,
    )

    best_params = boost_study.best_trial.params
    with open(conf.data.savedir + f"/{study_name}_best_params.pkl", "wb") as f:
        pickle.dump(best_params, f)

    console.print(f"Best parameters: {best_params}", justify="center")

    model = train_booster(conf, best_params, Xtrain, ytrain, Xtest, ytest)

    with open(conf.data.savedir + f"/{study_name}.pkl", "wb") as f:
        pickle.dump(model, f)

    console.log("Tuning Done", style="bold green", justify="center")

    console.save_html(conf.data.savedir + f"/{study_name}_run_report.html")


if __name__ == "__main__":
    run_command = " ".join([sys.executable, *sys.argv])
    console.log(f"[red]Run command[/red]: {run_command}")
    main()
