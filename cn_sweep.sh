#!/bin/bash
#SBATCH --partition=orion
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --job-name=cn-sweep
#SBATCH --output=cn_sweep_%j.out

# Clean stale files from $TMPDIR on all compute nodes (parallel via xargs).
# Removes: Rtmp*, tmp.*, nxf-*, nxf.*, ompi.*, singularity cache,
#           hsperfdata, libjansi, rstudio leftovers, resolv.conf,
#           build-temp-*, bundle-temp-* (Apptainer build cache)
# Usage: ./cn_sweep.sh   (or sbatch cn_sweep.sh)

mapfile -t NODES < <(sinfo -p orion -h -o "%n" | sort)

sweep_node() {
    local node=$1
    srun -n1 --partition=orion --nodelist="$node" --time=00:02:00 --cpus-per-task=1 --mem=512M bash -c '
        DIR="$(dirname "$TMPDIR")/$USER"

        if [ ! -d "$DIR" ]; then
            printf "%-8s  (no directory)\n" "$(hostname -s)"
            exit 0
        fi

        before=$(du -sh "$DIR" 2>/dev/null | cut -f1)
        count_before=$(find "$DIR" -maxdepth 1 -mindepth 1 2>/dev/null | wc -l)

        rm -rf "$DIR"/Rtmp* "$DIR"/tmp.* \
               "$DIR"/nxf-* "$DIR"/nxf.* \
               "$DIR"/ompi.* "$DIR"/singularity \
               "$DIR"/hsperfdata_* "$DIR"/libjansi* \
               "$DIR"/rstudio-server-* "$DIR"/resolv.conf \
               "$DIR"/build-temp-* "$DIR"/bundle-temp-* \
               "$DIR"/$USER 2>/dev/null

        after=$(du -sh "$DIR" 2>/dev/null | cut -f1)
        count_after=$(find "$DIR" -maxdepth 1 -mindepth 1 2>/dev/null | wc -l)

        printf "%-8s  %s (%d entries) -> %s (%d entries)\n" \
            "$(hostname -s)" "$before" "$count_before" "$after" "$count_after"
    ' 2>/dev/null
}
export -f sweep_node

printf '%s\n' "${NODES[@]}" | xargs -P "${#NODES[@]}" -I {} bash -c 'sweep_node "$@"' _ {}

# Author: Martin Paliocha <martpali> martin.paliocha@nmbu.no