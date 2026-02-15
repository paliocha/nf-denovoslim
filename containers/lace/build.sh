#!/bin/bash
# Build custom Lace container with .node -> .nodes patch

set -euo pipefail

CONTAINER_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$CONTAINER_DIR"

echo "Building Lace container with NetworkX .nodes patch..."
echo "Output: lace_1.14.1_patched.sif"

# Build with Apptainer (on a compute node to avoid overloading login node)
apptainer build lace_1.14.1_patched.sif Lace.def

echo "Done! Container: ${CONTAINER_DIR}/lace_1.14.1_patched.sif"
echo ""
echo "Update nextflow.config:"
echo "  lace_container = '\${projectDir}/containers/lace/lace_1.14.1_patched.sif'"
