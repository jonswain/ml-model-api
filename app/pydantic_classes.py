"""Pydantic classes for FastAPI validation and docs"""

from datetime import datetime
from pydantic import BaseModel, Field, conlist
from rdkit import Chem


class SmilesString(str):
    """Valid SMILES strings for Chemprop prediction"""

    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, smiles):
        """Checks submitted SMILES are valid else raises TypeError"""
        if not Chem.MolFromSmiles(smiles):
            raise TypeError(f"Invalid SMILES: {smiles}")
        return smiles


class RequestBody(BaseModel):
    """API request body format"""

    smiles: list[SmilesString] | None = Field(
        default=None,
        title="A list of valid SMILES strings",
        example=[
            "CC(=O)OC1=CC=CC=C1C(=O)O",
            "CN1C2CCC1C(C(C2)OC(=O)C3=CC=CC=C3)C(=O)OC",
        ],
    )


class Prediction(BaseModel):
    """Prediction format, uses processed Chemprop predictions"""

    smiles: str | None = Field(None, example="CC(=O)OC1=CC=CC=C1C(=O)O")
    prediction: float | None = Field(None, example=3.12)
    uncertainty: float | None = Field(None, example=0.58)


class Predictions(BaseModel):
    """API response format, a list of predictions"""

    predictions: conlist(Prediction)


class Details(BaseModel):
    """Detail format. Information about the machine learning model"""

    name: str | None = Field(None, example="LogD")
    version: str | None = Field(None, example="0.1.0")
    train_date: datetime | None = Field(
        None, example=datetime(year=2023, month=1, day=1)
    )
    model_type: str | None = Field(None, example="chemprop")
    data_source: str | None = Field(None, example="internal")
    model_reference: str | None = Field(None, example="trained internally")
    training_repo: str | None = Field(
        None, example="UPDATE"
    )
    prediction_type: str | None = Field(None, example="regression")
    features: str | None = Field(None, example="RDKit2D")
    calibration_method: str | None = Field(None, example="None")
    extra_information: str | None = Field(None, example="Extra information")
