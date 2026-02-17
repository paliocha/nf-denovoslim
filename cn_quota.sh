#!/bin/bash
#SBATCH --partition=orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --job-name=cn-quota
#SBATCH --output=cn_quota_%j.out

# Check /work (local TMPDIR) free space on all compute nodes.
# Usage: ./cn_quota.sh   (or sbatch cn_quota.sh)

mapfile -t NODES < <(sinfo -p orion -h -o "%n" | sort)

printf "%-8s %6s %6s %6s %5s  %-s\n" "NODE" "SIZE" "USED" "AVAIL" "USE%" "YOUR_USAGE"
printf "%s\n" "-------------------------------------------------------"

for node in "${NODES[@]}"; do
    srun -n1 --partition=orion --nodelist=$node --time=00:01:00 --cpus-per-task=1 --mem=512M bash -c '
        WORK_BASE="$(dirname "$TMPDIR")"
        read _ size used avail pct _ <<< $(df -h "$TMPDIR" 2>/dev/null | tail -1)
        my_usage=$(du -sh "$WORK_BASE/$USER" 2>/dev/null | cut -f1)
        [ -z "$my_usage" ] && my_usage="0"
        printf "%-8s %6s %6s %6s %5s  %s\n" "$(hostname -s)" "$size" "$used" "$avail" "$pct" "$my_usage"
    ' 2>/dev/null &
done
wait

# Author: Martin Paliocha <martpali> martin.paliocha@nmbu.no