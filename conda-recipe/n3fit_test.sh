#!/bin/bash
set -u
set -v
set -e

# Python tests for n3fit
pytest --pyargs n3fit
