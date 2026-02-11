#!/usr/bin/env python3
# filepath: /mnt/users/martpali/AnnualPerennial/nf-corset/bin/process_sample_sheet.py

import sys
import csv
import os
import re

def main():
    """
    Process a sample sheet CSV file to extract read pairs and generate grouping information.
    Format: sample,fastq_1,fastq_2,strandedness
    
    Creates two output files:
    1. read_pairs.txt: tab-delimited file with sample, R1, R2 paths
    2. grouping.txt: comma-separated list of group numbers for Corset
    """
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <sample_sheet.csv> <read_pairs_output> <grouping_output>")
        sys.exit(1)
    
    sample_sheet = sys.argv[1]
    read_pairs_output = sys.argv[2]
    grouping_output = sys.argv[3]
    
    samples = []
    
    # Read sample information from CSV
    with open(sample_sheet, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Extract sample name, fastq paths
            sample = row['sample'].strip()
            fastq_1 = row['fastq_1'].strip()
            fastq_2 = row['fastq_2'].strip()
            
            # Skip incomplete entries
            if not sample or not fastq_1 or not fastq_2:
                continue
                
            # Extract tissue and timepoint from sample name (format: BMED##_T#_[LR])
            match = re.match(r'([^_]+)_T([0-9])_([LR])', sample)
            if match:
                prefix = match.group(1)  # BMED##
                time_point = int(match.group(2))  # Number after T
                tissue = match.group(3)  # L or R
                
                # Store sample information
                samples.append({
                    'sample': sample,
                    'fastq_1': fastq_1,
                    'fastq_2': fastq_2,
                    'time_point': time_point,
                    'tissue': tissue
                })
    
    # Sort samples by timepoint then tissue
    samples.sort(key=lambda x: (x['time_point'], x['tissue']))
    
    # Generate grouping based on tissue and timepoint
    grouping = []
    
    for sample in samples:
        time_num = sample['time_point']
        if sample['tissue'] == 'L':
            # Leaf groups (T1_L=1, T2_L=3, T3_L=5, T4_L=7, T5_L=9)
            group = 2*time_num - 1
        else:
            # Root groups (T1_R=2, T2_R=4, T3_R=6, T4_R=8, T5_R=10)
            group = 2*time_num
        
        grouping.append(str(group))
    
    # Write the read pairs to output file
    with open(read_pairs_output, 'w') as f:
        for sample in samples:
            f.write(f"{sample['sample']}\t{sample['fastq_1']}\t{sample['fastq_2']}\n")
    
    # Write the grouping as a comma-separated list
    with open(grouping_output, 'w') as f:
        f.write(','.join(grouping))
    
    print(f"Processed {len(samples)} samples")
    print(f"Read pairs written to: {read_pairs_output}")
    print(f"Grouping written to: {grouping_output}")
    
    # Print groups for verification
    groups = {}
    for i, sample in enumerate(samples):
        group = grouping[i]
        if group not in groups:
            groups[group] = []
        groups[group].append(sample['sample'])
    
    print("\nGrouping summary:")
    for group in sorted(groups.keys(), key=int):
        print(f"Group {group}: {len(groups[group])} samples - {', '.join(groups[group])}")

if __name__ == "__main__":
    main()