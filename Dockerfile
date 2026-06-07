FROM mambaorg/micromamba:1.5.8

USER root

# System tools used by the pipeline / future external modules
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    curl \
    wget \
    git \
    ca-certificates \
    openjdk-17-jre-headless \
    hmmer \
    mafft \
    muscle \
    docker.io \
    tar \
    gzip \
    unzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy environment first for better Docker caching
COPY environment.yml /app/environment.yml

# Create the conda/mamba environment
RUN micromamba create -y -n compoundrank -f /app/environment.yml \
    && micromamba clean -a -y

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow

# Copy project code
COPY . /app

# Make sure the environment is active for commands
SHELL ["micromamba", "run", "-n", "compoundrank", "/bin/bash", "-c"]

# Install pip requirements too, in case environment.yml missed anything
RUN if [ -f requirements.txt ]; then pip install -r requirements.txt; fi

ENV PATH="/opt/conda/envs/compoundrank/bin:${PATH}"

CMD ["micromamba", "run", "-n", "compoundrank", "python", "-m", "uvicorn", "cpu_server.app:app", "--host", "0.0.0.0", "--port", "8000"]