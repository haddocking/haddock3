#==============================================================================================
FROM python:3.14-slim

# uv is used to speed up dependency resolution and installation
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

# mpi4py (mpi extra) needs an MPI implementation available at build time
RUN apt-get update \
  && apt-get install -y --no-install-recommends openmpi-bin libopenmpi-dev \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/software
COPY . .

RUN uv pip install --system --no-cache '.[mpi,dev,docs,notebooks]'

WORKDIR /data
ENTRYPOINT [ "haddock3" ]
#==============================================================================================
