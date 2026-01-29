from math import log
import cobra as cb
# import numpy as np
# from scipy.io import loadmat
from pathlib import Path
from cobra.sampling import OptGPSampler
import pandas as pd
import os 
import sys
import logging 


logging.basicConfig(format='%(asctime)s %(message)s')

logger = logging.getLogger("Sampling")
logger.setLevel(logging.INFO)
logger.addHandler((logging.StreamHandler(sys.stdout)))


path = Path("models/strain_models_gimme/")
savepath = Path("results/sampling_all/")

concordant_strains = pd.read_csv(
    "data/concordant_strains_ynb_gimme.csv",
)["Strain"]

concordant_strains = [strain.stem for strain in list(path.glob("*.mat"))]

if __name__ == "__main__":
    for strain in concordant_strains:
        # print(strain)
        if (savepath / f"{strain}.parquet").exists():
            
            continue
        logger.info(f"Sampling {strain}")
        model = cb.io.load_matlab_model(path / f"{strain}.mat")
        sampler = OptGPSampler(model, thinning=50, processes=4, seed=42)
        sample = sampler.sample(10000)
        sample.to_parquet(savepath / f"{strain}.parquet", compression="gzip")
        logger.info(f"{strain} done!")