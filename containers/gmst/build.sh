#!/usr/bin/env bash
# Build the GeneMarkS-T Apptainer container
#
# PREREQUISITES:
#   1. Download GeneMarkS-T from http://topaz.gatech.edu/GeneMark/license_download.cgi
#      - Select "GeneMarkS-T" and "Linux 64"
#      - Fill in the academic license form
#   2. Extract the tarball and copy gmst.pl into this directory:
#        tar xzf gmst_linux_64.tar.gz
#        cp gmst_linux_64/gmst.pl containers/gmst/
#   3. Run this script from the containers/gmst/ directory
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DEF="${SCRIPT_DIR}/GeneMarkST.def"
SIF="${SCRIPT_DIR}/gmst.sif"

if [[ ! -f "${SCRIPT_DIR}/gmst.pl" ]]; then
    echo "ERROR: gmst.pl not found in ${SCRIPT_DIR}/" >&2
    echo "" >&2
    echo "Download GeneMarkS-T from http://topaz.gatech.edu/GeneMark/license_download.cgi" >&2
    echo "Extract and copy gmst.pl here before building." >&2
    exit 1
fi

echo "Building ${SIF} from ${DEF} ..."
apptainer build --fakeroot "${SIF}" "${DEF}"
echo "Done. Test with:  apptainer exec ${SIF} gmst.pl 2>&1 | head -5"
