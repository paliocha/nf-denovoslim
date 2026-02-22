#!/usr/bin/env bash
# Build the TransAnnot Apptainer container (eggNOG v7 branch)
# Run from this directory on an Orion login node
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DEF="${SCRIPT_DIR}/Transannot.def"
SIF="${SCRIPT_DIR}/transannot_4.0.0_eggnog7.sif"

echo "Building ${SIF} from ${DEF} ..."
apptainer build --fakeroot "${SIF}" "${DEF}"
echo "Done. Test with:  apptainer exec ${SIF} transannot -h"
