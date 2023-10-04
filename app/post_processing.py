"""Functions to process chemprop predictions"""
import numpy as np


def process_predictions(results_df):
    """Processes chemprop regression predictions

    Args:
        results_df (DataFrame): A DataFrame containing a column of raw
        chemprop classification predictions

    Returns:
        DataFrame: Returns predictions, will add uncertainty at a later stage
    """
    return results_df
    