"""FastAPI for making predictions using a machine learning model"""

import os
from fastapi import FastAPI, Request, APIRouter
from fastapi.responses import RedirectResponse
from app.pydantic_classes import RequestBody, Predictions, Details
from app.prediction import make_prediction
from . import config

app = FastAPI(
    root_path=os.getenv("DEPLOY_PATH", ""),
    title=f"{config.MODEL_DETAIL['name']} prediction model",
    description=config.MODEL_DESCRIPTION,
    version=config.MODEL_DETAIL["version"],
)
router = APIRouter()
app.include_router(router)


@app.get("/", include_in_schema=False)
async def docs_redirect():
    """Redirects root to docs"""
    return RedirectResponse(url="/docs")


@app.get("/healthcheck", tags=["API endpoints"])
def healthcheck(request: Request) -> dict:
    """Healthcheck API endpoint, returns message if API is running"""
    return {
        "message": "API is running",
        "root_path": request.scope.get("root_path"),
    }


@app.post(
    "/predict",
    response_model=Predictions,
    response_model_exclude_none=True,
    tags=["API endpoints"],
)
def predict(data: RequestBody) -> Predictions:
    """Prediction API endpoint, receives a list of SMILES strings and returns predictions"""
    predictions = make_prediction(smiles=data.smiles)
    return {"predictions": predictions}


@app.get("/detail", tags=["API endpoints"])
def detail() -> Details:
    """Detail API endpoint, returns details of model as a dictionary"""
    return config.MODEL_DETAIL
