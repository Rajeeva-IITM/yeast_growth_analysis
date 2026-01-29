import pickle
import lightgbm as lgb
from omegaconf import DictConfig
import hydra
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from utils import get_data, get_model_paths
from rich.console import Console
from dotenv import load_dotenv
from pathlib import Path

load_dotenv()


def train_booster(config: DictConfig, params, Xtrain, ytrain, Xtest, ytest):
    constant_params = {  # These parameters won't change regardless of the trial and model
        "objective": config.model_params.objective,
        "verbosity": 1,
        "force_col_wise": config.model_params.force_col_wise,
        "early_stopping_rounds": config.model_params.early_stopping_rounds,
        "num_threads": config.model_params.num_threads,
        "boosting_type": config.model_params.boosting_type,
        "device_type": config.model_params.device_type,
        "gpu_use_dp": config.model_params.gpu_use_dp,
        "seed": config.seed,
    }

    final_params = constant_params | params

    model = lgb.train(
        params=final_params,
        train_set=lgb.Dataset(Xtrain, label=ytrain),
        valid_sets=[lgb.Dataset(Xtest, label=ytest)],
    )

    return model


@hydra.main(config_path="../configs/", version_base="1.3", config_name="conf")
def main(conf: DictConfig):
    console = Console()

    study_name = (
        conf.models.Boosting.name
        if conf.model_params.boosting_type == "gbdt"
        else conf.models.RandomForest.name
    )

    # locating params
    model_param_paths = get_model_paths(
        run_path=conf.model_load_keys.run_path,
        model_type=conf.model_load_keys.model_type + "_best_params",
        prefix=conf.model_load_keys.prefix,
        suffix=conf.model_load_keys.suffix,
        model_names=conf.model_load_keys.model_names,
    )
    console.print(model_param_paths)

    # loading parameters
    for name, paths in model_param_paths.items():
        # Data processing
        console.log("Processing data", style="bold red", justify="center")

        Xtrain, ytrain = get_data(
            conf.data_paths.get(name), run_type=conf.run_type, return_as_Xy=True
        )
        Xtrain, Xtest, ytrain, ytest = train_test_split(
            Xtrain, ytrain, test_size=conf.testing.test_frac, random_state=conf.seed
        )

        scaler = StandardScaler()
        Xtrain = scaler.fit_transform(Xtrain)
        Xtest = scaler.transform(Xtest)
        console.log("Data processed", style="bold green", justify="center")

        for model_param_path in paths:
            model_param_path = Path(model_param_path)
            console.print(model_param_path)

            with open(model_param_path, "rb") as f:
                params = pickle.load(f)

            model = train_booster(
                config=conf,
                params=params,
                Xtrain=Xtrain,
                ytrain=ytrain,
                Xtest=Xtest,
                ytest=ytest,
            )

            with open(model_param_path.parent / f"{study_name}.pkl", "wb") as f:
                pickle.dump(model, f)

    return None


if __name__ == "__main__":
    main()
