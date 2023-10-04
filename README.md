# LogS-model-API

A FastAPI endpoint for a machine learning model to predict LogS.
  
## Local development

The following commands will setup an environment where you can run and test the application locally:

```shell
git clone UPDATE
cd logs-model-api
mamba env create -f envs/dev.yml
conda activate LogS-model-API
code .
```

### Running a local uvicorn ASGI web server

A local server, running on the default port (8000), can be started as follows:

```shell
uvicorn app.main:app
```

Once started, you can browse to <http://localhost:8000> to test the app.

### Running a local gunicorn WSGI web server in docker

A docker container running gunicorn can be started locally using docker compose:

```shell
docker compose up -d
```

The base image handles sensible gunicorn configuration options.

The `backend` service can be accessed at port 5000 (<http://localhost:5000>)

## Directory structure

    ├── README.md               <- The top-level README for developers using this project.
    ├── app
    │   ├── __init__.py         <- Makes LogS-model-API a Python module.
    │   ├── config.py           <- Config data for machine learning model API.
    │   ├── main.py             <- API functions for hosting.
    │   ├── post_processing.py  <- Scripts for processing predictions.
    │   ├── prediction.py       <- Scripts for making predictions.
    │   └── pydantic_classes.py <- Pydantic classes for API data validation.
    │
    ├── envs                    <- Conda environments.
    │   └── dev.yml             <- Conda environment for development.
    │
    ├── model                   <- Directory to place your machine learning model in.
    │
    ├── tests                   <- Test functions for deployment.
    │           
    ├── .dockerignore
    │           
    ├── .gitignore
    │           
    ├── Dockerfile              <- Docker compose for development.
    │
    └── requirements.txt        <- Requirements file for dependencies in Docker container.
