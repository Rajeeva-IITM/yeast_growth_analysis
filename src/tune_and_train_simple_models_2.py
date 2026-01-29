import logging
import pickle

# from pathlib import Path
from typing import Callable

import hydra
import numpy as np
import polars as pl
import optuna
from dotenv import load_dotenv
from sklearn.discriminant_analysis import (
    LinearDiscriminantAnalysis,
    QuadraticDiscriminantAnalysis,
)
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.model_selection import KFold, train_test_split
from sklearn.preprocessing import StandardScaler
from omegaconf import DictConfig
from rich.console import Console
from rich.logging import RichHandler

load_dotenv()

console = Console(record=True)

logger = logging.getLogger()
logger.addHandler(RichHandler(console=console))

optuna.logging.disable_default_handler()
optuna.logging.enable_propagation()


def tune_LDA(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the Linear Discriminant Analysis model."""
    solver = trial.suggest_categorical("solver", ["svd", "lsqr", "eigen"])
    shrinkage = (
        trial.suggest_float("shrinkage", 0.0, 1.0)
        if solver in ["lsqr", "eigen"]
        else None
    )

    params = {"solver": solver, "shrinkage": shrinkage}

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )
    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = LinearDiscriminantAnalysis(**params)
        model.fit(X_train, y_train)
        scores.append(
            hydra.utils.call(config.metric, _args_=(y_test, model.predict(X_test)))
        )

    return np.mean(scores)


def tune_QDA(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the Quadratic Discriminant Analysis model."""
    reg_param = trial.suggest_float("reg_param", 0.0, 1.0)

    params = {"reg_param": reg_param}

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )
    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = QuadraticDiscriminantAnalysis(**params)
        model.fit(X_train, y_train)
        scores.append(
            hydra.utils.call(config.metric, _args_=(y_test, model.predict(X_test)))
        )

    return np.mean(scores)


def tune_NearestNeighbors(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the K-Nearest Neighbors classifier."""
    n_neighbors = trial.suggest_int("n_neighbors", 1, 50)
    weights = trial.suggest_categorical("weights", ["uniform", "distance"])
    # metric = trial.suggest_categorical("metric", ["euclidean", "manhattan", "minkowski"])

    params = {"n_neighbors": n_neighbors, "weights": weights}

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )
    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = KNeighborsClassifier(**params)
        model.fit(X_train, y_train)
        scores.append(
            hydra.utils.call(config.metric, _args_=(y_test, model.predict(X_test)))
        )

    return np.mean(scores)


def tune_NaiveBayes(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the Gaussian Naive Bayes classifier."""
    var_smoothing = trial.suggest_float("var_smoothing", 1e-9, 1e-7, log=True)

    params = {"var_smoothing": var_smoothing}

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )
    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = GaussianNB(**params)
        model.fit(X_train, y_train)
        scores.append(
            hydra.utils.call(config.metric, _args_=(y_test, model.predict(X_test)))
        )

    return np.mean(scores)


def tune_GaussianProcesses(trial: optuna.Trial, X, y, config: DictConfig):
    """Tunes the Gaussian Process Classifier."""
    kernel = RBF(length_scale=trial.suggest_float("length_scale", 1e-2, 1e2, log=True))

    params = {"kernel": kernel}

    splits = KFold(
        n_splits=config.kfold_params.n_splits,
        shuffle=config.kfold_params.shuffle,
        random_state=config.kfold_params.seed,
    )
    scores = []
    for train_idx, test_idx in splits.split(X, y):
        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]
        model = GaussianProcessClassifier(**params)
        model.fit(X_train, y_train)
        scores.append(
            hydra.utils.call(config.metric, _args_=(y_test, model.predict(X_test)))
        )

    return np.mean(scores)


def run_study(
    study_name: str, tune_objective: Callable, conf: DictConfig, X, y, n_jobs: int = 1
):
    """Runs an Optuna study to optimize hyperparameters for a given model."""
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
        catch=(ValueError,),
        n_jobs=n_jobs,
    )

    return study


@hydra.main(version_base=None, config_path="../configs/", config_name="conf")
def main(conf: DictConfig):
    """Main function to process data and tune models."""
    console.log("Processing data", style="bold red", justify="center")

    df = pl.read_ipc(conf.data.path)

    X = (
        df.drop(["Phenotype", "Condition", "Strain"])
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

    lda_study = run_study(
        study_name=conf.models.LDA.name,
        tune_objective=tune_LDA,
        conf=conf,
        X=Xtrain,
        y=ytrain,
    )
    best_model = LinearDiscriminantAnalysis(**lda_study.best_trial.params)
    best_model.fit(Xtrain, ytrain)
    pickle.dump(best_model, open(conf.models.LDA.model_savename, "wb"))

    console.log("LDA Done", style="bold green", justify="center")

    # QDA
    qda_study = run_study(
        study_name=conf.models.QDA.name,
        tune_objective=tune_QDA,
        conf=conf,
        X=Xtrain,
        y=ytrain,
    )
    best_model = QuadraticDiscriminantAnalysis(**qda_study.best_trial.params)
    best_model.fit(Xtrain, ytrain)
    pickle.dump(best_model, open(conf.models.QDA.model_savename, "wb"))

    console.log("QDA Done", style="bold green", justify="center")

    # NB
    nb_study = run_study(
        study_name=conf.models.NB.name,
        tune_objective=tune_NaiveBayes,
        conf=conf,
        X=Xtrain,
        y=ytrain,
    )
    best_model = GaussianNB(**nb_study.best_trial.params)
    best_model.fit(Xtrain, ytrain)
    pickle.dump(best_model, open(conf.models.NB.model_savename, "wb"))

    console.log("NB Done", style="bold green", justify="center")

    # NearestNeighbors

    knn_study = run_study(
        study_name=conf.models.KNN.name,
        tune_objective=tune_NearestNeighbors,
        conf=conf,
        X=Xtrain,
        y=ytrain,
    )
    best_model = KNeighborsClassifier(**knn_study.best_trial.params)
    best_model.fit(Xtrain, ytrain)
    pickle.dump(best_model, open(conf.models.KNN.model_savename, "wb"))

    console.log("KNN Done", style="bold green", justify="center")

    # # GP
    # gp_study = run_study(
    #     study_name=conf.models.GP.name,
    #     tune_objective=tune_GaussianProcesses,
    #     conf=conf,
    #     X=Xtrain,
    #     y=ytrain,
    # )
    # best_model = GaussianProcessClassifier(**gp_study.best_trial.params)
    # best_model.fit(Xtrain, ytrain)
    # pickle.dump(best_model, open(conf.models.GP.model_savename, "wb"))

    # console.log("GP Done", style="bold green", justify="center")


if __name__ == "__main__":
    main()
