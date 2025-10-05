#!/usr/bin/env bash

set -euo pipefail

# directory containing fastq read   
input_path="../data/trimmed_reads/"

# directory containing aligned reads (bam files)
output_path="../data/mapping/bam_unsorted/"

# directory containing mapping statistics for each sample
stats_path="../data/mapping/bam_stats/"

# reference genome
ref_genome="../../../refdata/hg38/hg38.fasta

# make output and mapping statistics directories if they don't exist
mkdir -p "${output_path}" "${stats_path}"

# checks if genome index exists and creates index if it doesn't exist
if [[ ! -f "${ref_genome}.amb" || \
      ! -f "${ref_genome}.ann" || \
      ! -f "${ref_genome}.bwt" || \
      ! -f "${ref_genome}.pac" || \
      ! -f "${ref_genome}.sa" ]]
then
    echo "Creating genome index for reference genome..."
    bwa index "${ref_genome}"
fi

# loops through the trimmed fastq files
for file in "${input_path}"*_trimmed_1.fastq.gz
do
    # remove the suffix (_trimmed_1.fastq.gz) from file name
    # example (sample_trimmed_1.fastq.gz becomes sample)
    sample=$(basename "${file}" _trimmed_1.fastq.gz)

    echo "Processing ${sample}..."

    # Assigns Read Group
    RG="@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA"

    # run BWA and output BAM
    bwa mem \
    -t 8 \
    -R "${RG}" \
    "${ref_genome}" \
    "${input_path}${sample}_trimmed_1.fastq.gz" \
    "${input_path}${sample}_trimmed_2.fastq.gz" \
    | samtools view -b -o "${output_path}${sample}.bam"

    # run samtools flagstat
    samtools flagstat "${output_path}${sample}.bam" > "${stats_path}${sample}.flagstat.txt"

    echo "Completed ${sample}"
done

echo "All mapping done"

