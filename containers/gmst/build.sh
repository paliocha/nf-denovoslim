#!/usr/bin/env bash
# Build the GeneMarkS-T Apptainer container
#
# PREREQUISITES:
#   1. Download GeneMarkS-T from http://topaz.gatech.edu/GeneMark/license_download.cgi
#      - Select "GeneMarkS-T" and "Linux 64"
#      - Also download the 64-bit license key (gm_key_64.gz)
#      - Fill in the academic license form
#   2. Extract files in this directory:
#        tar xzf gmst_linux_64.tar.gz
#        gunzip -k gm_key_64.gz
#   3. Run this script from the containers/gmst/ directory
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DEF="${SCRIPT_DIR}/GeneMarkST.def"
SIF="${SCRIPT_DIR}/gmst.sif"

REQUIRED_FILES=(gmst.pl gmhmmp probuild Gibbs3 gm_key_64)
MISSING=0
for f in "${REQUIRED_FILES[@]}"; do
    if [[ ! -f "${SCRIPT_DIR}/${f}" ]]; then
        echo "ERROR: ${f} not found in ${SCRIPT_DIR}/" >&2
        MISSING=1
    fi
done

if [[ $MISSING -eq 1 ]]; then
    echo "" >&2
    echo "Download GeneMarkS-T from http://topaz.gatech.edu/GeneMark/license_download.cgi" >&2
    echo "Extract gmst_linux_64.tar.gz and gunzip gm_key_64.gz in this directory." >&2
    exit 1
fi

echo "Building ${SIF} from ${DEF} ..."
apptainer build --fakeroot "${SIF}" "${DEF}"
echo "Done. Test with:  apptainer exec ${SIF} gmst.pl 2>&1 | head -5"
