FROM python:3.10 AS app

WORKDIR /app
ENV PATH="/app/.venv/bin:$PATH"
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/
COPY --chown=1000:100 main /app/main/
COPY --chown=1000:100 pyproject.toml /app/pyproject.toml

RUN uv sync && uv pip install jupyterlab
EXPOSE 8888
VOLUME "/app/main/data/local_files"
CMD ["jupyter", "lab", "--allow-root", "--no-browser", "--ip=0.0.0.0", "--port=8888", "--notebook-dir=/app/main", "--NotebookApp.token=''", "--NotebookApp.password=''"]
