"""Functions to make predictions using machine a learning model"""

import pickle
import pandas as pd
from molfeat.calc import FPCalculator, get_calculator
from molfeat.trans import MoleculeTransformer
from rdkit import Chem
from .post_processing import process_predictions


def generate_features(smiles):
    """Generates Morgan count fingerprints and all RDKit 2D descriptors.

    Args:
        smiles (list): List of canonical smiles.

    Returns:
        Pandas dataframe: The combined features dataframe.
    """

    params = {
        "radius": 3,
        "nBits": 2048,
        "useChirality": True,
        "useFeatures": True,
    }

    fp_calc = FPCalculator("ecfp-count", **params)
    fp_transf = MoleculeTransformer(fp_calc, n_jobs=-1)

    rdkit_calc = get_calculator("desc2d")
    rdkit_transf = MoleculeTransformer(rdkit_calc, n_jobs=-1)

    df_desc = pd.DataFrame(fp_transf(smiles), columns=fp_calc.columns)
    df_rdkit = pd.DataFrame(rdkit_transf(smiles), columns=rdkit_calc.columns)
    df_rdkit.drop("Alerts", axis=1, inplace=True)

    return pd.concat([df_desc, df_rdkit], axis=1)


def make_prediction(smiles):
    """Makes predictions using sklearn model and SMILES provided

    Args:
        smiles (List): A list of valid SMILES strings

    Returns:
        List[Dict]: A list of the processed predictions of each given SMILES string
    """
    # Prepare SMILES
    results_df = pd.DataFrame(smiles, columns=["smiles"])
    results_df["canonical_smiles"] = results_df.smiles.apply(Chem.CanonSmiles)
    canonical_smiles = set(results_df.canonical_smiles.to_list())
    canonical_df = pd.DataFrame(canonical_smiles, columns=["canonical_smiles"])

    # Load model:
    with open("model/model.pkl", "rb") as file:
        regressor = pickle.load(file)

    # Make prediction
    preds = regressor.predict(generate_features(canonical_smiles))
    canonical_df["prediction"] = preds

    # Process predictions
    processed_df = process_predictions(canonical_df)

    results_df = results_df.merge(processed_df, how="left", on="canonical_smiles")

    return results_df.to_dict("records")
