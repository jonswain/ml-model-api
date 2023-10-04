"""Config file for LogS model API"""

from datetime import datetime

MODEL_DESCRIPTION = """A machine learning model to predict LogS.
* Data source: AI4SD Summer School
* Training repo: https://github.com/jonswain/solubility_prediction"""

MODEL_DETAIL = {
    "name": "LogS",
    "version": "0.1.0",
    "train_date": datetime(year=2023, month=7, day=13),
    "model_type": "sklearn",
    "data_source": "AI4SD Summer School",
    "training_repo:": "https://github.com/jonswain/solubility_prediction",
    "prediction_type": "regression",
    "features": "ECFP3 and RDKit2D",
    "extra_information": None,
}
