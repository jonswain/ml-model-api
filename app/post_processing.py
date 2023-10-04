"""Functions to process predictions"""

import pickle
import heapq
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

# Initialize objects
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

with open("data/fps.pkl", "rb") as file:
    train_fps = pickle.load(file)


def mean_highest_tanimoto(smiles):
    """Calculates a confidence value for a prediction
    Confidence is calculated as the mean of the 5 highest tanimoto similarities
    from the training set.

    Args:
        smiles (str): SMILES string for prediction

    Returns:
        float: Mean of the 5 highest tanimoto similarities in the training set
    """
    mol = Chem.MolFromSmiles(smiles)
    fp = fpgen.GetFingerprint(mol)
    mean_similarity = (
        sum(heapq.nlargest(5, DataStructs.BulkTanimotoSimilarity(fp, train_fps))) / 5
    )

    return mean_similarity


def process_predictions(results_df):
    """Processes predictions, adds a confidence value for each prediction
    Confidence is calculated as the mean of the 5 highest tanimoto similarities
    from the training set.

    Args:
        results_df (DataFrame): A DataFrame containing a column of raw
        chemprop classification predictions

    Returns:
        DataFrame: Returns predictions, with added confidence column
    """
    results_df["confidence"] = results_df.canonical_smiles.apply(mean_highest_tanimoto)

    return results_df
