FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PATH=/opt/gnina:${PATH}

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    openbabel \
    fpocket \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt pyproject.toml /app/
COPY compoundrank /app/compoundrank
COPY run_pipeline.py /app/

RUN python3 -m pip install --break-system-packages --no-cache-dir .

# GNINA is GPU/CUDA-specific and is intentionally not baked into this image.
# Mount a working GNINA binary at /opt/gnina/gnina or extend this image.
ENTRYPOINT ["compoundrank"]
