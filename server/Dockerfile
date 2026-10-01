# Lightweight Google Cloud Run Container for HADDOCK3 FastAPI Control Plane
FROM python:3.11-slim

# Prevent Python from buffering stdout/stderr and writing bytecode
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PORT=8080

WORKDIR /app

# Install minimal OS dependencies for network operations
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install FastAPI server dependencies and Modal client
RUN pip install --no-cache-dir \
    fastapi==0.142.2 \
    "uvicorn[standard]>=0.30.0" \
    pydantic==2.13.0 \
    pydantic-settings>=2.0.0 \
    python-multipart>=0.0.9 \
    httpx>=0.28.0 \
    modal>=1.3.4

# Copy server code
COPY server/ /app/server/

# Ensure workspace directory exists
RUN mkdir -p /app/server/workspace

EXPOSE 8080

# Launch Uvicorn bound to Cloud Run PORT environment variable
CMD exec uvicorn server.main:app --host 0.0.0.0 --port ${PORT:-8080}
