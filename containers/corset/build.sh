#!/usr/bin/env bash
# Build the Corset Apptainer container
# Run from this directory on an Orion login node
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DEF="${SCRIPT_DIR}/Corset.def"
SIF="${SCRIPT_DIR}/corset_1.10.sif"

echo "Building ${SIF} from ${DEF} ..."
apptainer build --fakeroot "${SIF}" "${DEF}"
echo "Done. Test with:  apptainer exec ${SIF} corset -t 4"
