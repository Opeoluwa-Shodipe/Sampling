#!/usr/bin/env bash

set -euo pipefail

# directory containing fastq read   
input_path="../data/raw_reads/"

# directory containing trimmed fastq reads
output_path="../data/trimmed_reads/"

# make output directory if it doesn't exist
mkdir -p "${output_path}"

for file in "${input_path}"*_1.fastq.gz
do
    # remove the suffix (_1.fastq.gz) from file name
    # example (sample_1.fastq.gz becomes sample)
    sample=$(basename "${file}" _1.fastq.gz)

    echo "Processing ${sample}..."

    # run fastp
    fastp \
    -i "${input_path}${sample}"_1.fastq.gz \
    -I "${input_path}${sample}"_2.fastq.gz \
    -o "${output_path}${sample}"_trimmed_1.fastq.gz \
    -O "${output_path}${sample}"_trimmed_2.fastq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 30 \
    --length_required  50 \
    --html "${output_path}${sample}"_fastp_report.html \
    --json "${output_path}${sample}"_fastp_report.json 

    echo "Completed ${sample}"
done

echo "Trimming completed successfully"
