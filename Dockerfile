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
    tar \
    gzip \
    unzip \
    procps \
    xz-utils \
    && rm -rf /var/lib/apt/lists/*

# Install Docker CLI only.
# We do NOT need the Docker daemon inside this container.
# The host Docker socket is mounted at runtime.
RUN curl -L https://download.docker.com/linux/static/stable/x86_64/docker-27.5.1.tgz -o /tmp/docker.tgz \
    && tar -xzf /tmp/docker.tgz -C /tmp \
    && mv /tmp/docker/docker /usr/local/bin/docker \
    && chmod +x /usr/local/bin/docker \
    && rm -rf /tmp/docker /tmp/docker.tgz

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
ENV PATH="/usr/local/bin:/usr/bin:/bin:${PATH}"

EXPOSE 8000

CMD ["python", "-m", "uvicorn", "cpu_server.app:app", "--host", "0.0.0.0", "--port", "8000"]