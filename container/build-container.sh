#!/bin/bash

set -euo pipefail

pkgs=(
    build-essential
    ca-certificates
    curl
    gfortran
    git
)

apt-get update -y
apt-get install -y --no-install-recommends "${pkgs[@]}"
apt-get clean
rm -rf /var/lib/apt/lists/*

# install micromamba
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -C /usr/local/bin --strip-components=1 -xvj bin/micromamba

# install Python, LHAPDF, and pandoc into the base environment
micromamba install --yes \
    --root-prefix="${MAMBA_ROOT_PREFIX}" \
    --name base \
    --channel conda-forge \
    "python=${PYTHON_V}" \
    lhapdf \
    pandoc \
    pip

# clean package cache
micromamba clean --all --yes
