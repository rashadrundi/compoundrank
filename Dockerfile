FROM python:3.12-slim

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN python -m pip install --upgrade pip && \
    python -m pip install -r /app/requirements.txt

COPY homolog_search.py /app/homolog_search.py
COPY folding.py /app/folding.py
COPY run_pipeline.py /app/run_pipeline.py

CMD ["python", "run_pipeline.py"]