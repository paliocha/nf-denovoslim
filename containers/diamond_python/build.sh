#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
apptainer build diamond_python.sif Diamond_Python.def
echo "Built: diamond_python.sif"
