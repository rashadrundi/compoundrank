FROM python:3.12-slim

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    curl \
    wget \
    git \
    ca-certificates \
    default-jre-headless \
    hmmer \
    mafft \
    muscle \
    docker.io \
    tar \
    gzip \
    unzip \
    procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/requirements.txt

RUN python -m pip install --upgrade pip \
    && python -m pip install -r /app/requirements.txt

RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow

COPY . /app

ENV COMPOUNDRANK_REPO_ROOT=/app
ENV COMPOUNDRANK_DEPLOY_ROOT=/opt/compoundrank
ENV COMPOUNDRANK_DATA_ROOT=/opt/compoundrank/data
ENV COMPOUNDRANK_JOBS_ROOT=/opt/compoundrank/jobs
ENV INTERPRO_DATA_ROOT=/opt/compoundrank/data/interpro_data
ENV VOGDB_DATA_ROOT=/opt/compoundrank/data/vogdb_data

EXPOSE 8000

CMD ["python", "-m", "uvicorn", "cpu_server.app:app", "--host", "0.0.0.0", "--port", "8000"]