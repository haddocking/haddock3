#!/usr/bin/env bash
set -e

# Run the HADDOCK3 FastAPI Cloud Server locally
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

export PYTHONPATH="$REPO_ROOT:$PYTHONPATH"
PORT="${PORT:-8080}"
HOST="${HOST:-0.0.0.0}"

echo "=========================================================="
echo " Starting HADDOCK3 GPU Serverless Control Plane"
echo " Host: http://$HOST:$PORT"
echo " Docs: http://$HOST:$PORT/docs"
echo "=========================================================="

if [ -f "$REPO_ROOT/.venv/bin/uvicorn" ]; then
    exec "$REPO_ROOT/.venv/bin/uvicorn" server.main:app --host "$HOST" --port "$PORT" --reload
else
    exec uvicorn server.main:app --host "$HOST" --port "$PORT" --reload
fi
