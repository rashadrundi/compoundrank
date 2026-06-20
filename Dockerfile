FROM nvidia/cuda:12.8.1-cudnn-runtime-ubuntu24.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MINIFORGE_VERSION=25.9.1-0
ARG GNINA_VERSION=1.3.2

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    CONDA_DIR=/opt/conda \
    PATH=/opt/conda/bin:/usr/local/bin:${PATH} \
    XDG_CACHE_HOME=/cache \
    COLABFOLD_CACHE=/cache/colabfold \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Core operating-system dependencies.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    wget \
    git \
    bzip2 \
    openbabel \
    tini \
    libgomp1 \
    libglib2.0-0 \
    libgl1 \
    libxext6 \
    libxrender1 \
    libsm6 \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge and the native ColabFold utilities.
RUN curl -fsSL \
    "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh" \
    -o /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p "${CONDA_DIR}" \
    && rm /tmp/miniforge.sh \
    && conda install -y \
        -c conda-forge \
        -c bioconda \
        python=3.13 \
        pip \
        kalign2=2.04 \
        hhsuite=3.3.0 \
        mmseqs2=18.8cc5c \
    && conda clean -afy

# Install GPU-enabled ColabFold.
RUN python -m pip install --upgrade pip setuptools wheel \
    && python -m pip install \
        "colabfold[alphafold,openmm]==1.6.1" \
        "jax[cuda]<0.8" \
        "openmm[cuda12]"

# Install the CUDA 12.8 GNINA binary.
RUN curl -fL \
    "https://github.com/gnina/gnina/releases/download/v${GNINA_VERSION}/gnina.${GNINA_VERSION}.cuda12.8" \
    -o /usr/local/bin/gnina \
    && chmod 0755 /usr/local/bin/gnina

# Persistent working locations.
RUN mkdir -p \
    /cache/colabfold \
    /work/input \
    /work/jobs \
    /work/output

WORKDIR /app

# Copy the complete pulled repository so configuration files and package
# resources are not accidentally omitted.
COPY . /app

# Install CompoundRank/EXORCIST from pyproject.toml.
RUN python -m pip install . \
    && python -m pip check \
    && command -v compoundrank \
    && command -v colabfold_batch \
    && command -v gnina

ENTRYPOINT ["/usr/bin/tini", "--", "compoundrank"]