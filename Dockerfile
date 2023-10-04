FROM tiangolo/uvicorn-gunicorn-fastapi:python3.11-slim

WORKDIR /app

COPY requirements.txt requirements.txt

# Required to install descriptastorus
RUN apt-get -y update
RUN apt-get -y install git

# Fix ImportError
RUN apt-get install libxrender1

RUN pip install --no-cache-dir --upgrade -r requirements.txt

COPY app app
COPY tests tests
COPY model model

RUN apt-get update && apt-get install -y \
  cifs-utils keyutils \
  && rm -rf /var/lib/apt/lists/*
