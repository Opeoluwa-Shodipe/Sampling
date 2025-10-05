#!/usr/bin/env bash

set -euo pipefail

# directory containing raw, unsorted bam files   
bam_path="../data/mapping/bam_unsorted/"

# reference genome
ref_genome="../../../refdata/hg38/hg38.fasta"

# directory to store sorted bam files
sorted_bam_path="../data/mapping/bam_sorted/"

# directory to store marked duplicate BAM files
marked_bam_path="../data/mapping/bam_marked/"   

# directory to store recalibration table
recal_table_path="../data/mapping/recal_tables/"

# directory to store recalibrated bam files
recal_bam_path="../data/mapping/bam_recal/" 

# Path to known variant sites (dbSNP) used in base recalibration
dbsnp_path="../../../refdata/hg38/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

# Path to known indel sites used in base recalibration
known_indels_path="../../../refdata/hg38/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"

# Base variant directory
variant_base="../data/variants/"

# Directory to store raw per-sample GVCF files (input for joint genotyping)
per_sample_gvcf_path="${variant_base}per_sample_gvcfs/"


# create directories if they don't exist
mkdir -p \
    "${sorted_bam_path}" "${marked_bam_path}" "${recal_table_path}" \
    "${recal_bam_path}" "${variant_base}" "${per_sample_gvcf_path}"  

# Ensure reference genome index exist 
if [[ ! -f "${ref_genome}.fai" ]]; then
    echo "Creating FASTA index for reference genome..."
    samtools faidx "${ref_genome}"
fi

# Creates sequence dictionary if it doesn't exist and skip the creation of the sequence dictionary if it does exist
if [[ ! -f "${ref_genome}.dict" ]]; then
    echo "Creating sequence dictionary for reference genome..."
    gatk CreateSequenceDictionary \
        -R "${ref_genome}" \
        -O "${ref_genome}.dict"
fi

# Loop through the bam files
for file in "${bam_path}"*.bam
do
    # strip out the .bam in the file
    sample="$(basename "${file}" .bam)"

    echo "Processing sample: ${sample}"

    # sort bam files
    echo "Sorting BAM for ${sample}..."
    gatk SortSam \
        -I "${bam_path}${sample}.bam" \
        -O "${sorted_bam_path}${sample}_sorted.bam" \
        -SORT_ORDER coordinate

    # index the sorted BAM file
    echo "Indexing sorted BAM for ${sample}..."
    gatk BuildBamIndex \
        -I "${sorted_bam_path}${sample}_sorted.bam" \
        -O "${sorted_bam_path}${sample}_sorted.bai"

    # mark duplicates
    echo "Marking duplicates for ${sample}..."
    gatk MarkDuplicates \
        -I "${sorted_bam_path}${sample}_sorted.bam" \
        -O "${marked_bam_path}${sample}_marked.bam" \
        -M "${marked_bam_path}${sample}_metrics.txt"
        
    # index the marked BAM
    echo "Indexing marked BAM for ${sample}..."
    gatk BuildBamIndex \
        -I "${marked_bam_path}${sample}_marked.bam" \
        -O "${marked_bam_path}${sample}_marked.bai"

    # Base Quality Score Recalibration (build recalibration table)
    echo "Running BaseRecalibrator for ${sample}..."
    gatk BaseRecalibrator \
        -I "${marked_bam_path}${sample}_marked.bam" \
        -R "${ref_genome}" \
        --known-sites "${dbsnp_path}" \
        --known-sites "${known_indels_path}" \
        -O "${recal_table_path}${sample}_recal.table"

    # Apply Base Quality Score Recalibration
    echo "Applying BQSR for ${sample}..."
    gatk ApplyBQSR \
        -R "${ref_genome}" \
        -I "${marked_bam_path}${sample}_marked.bam" \
        --bqsr-recal-file "${recal_table_path}${sample}_recal.table" \
        -O "${recal_bam_path}${sample}_recal.bam"

    # Call Variants
    echo "Calling variants for ${sample}..."
    gatk HaplotypeCaller \
        -I "${recal_bam_path}${sample}_recal.bam" \
        -R "${ref_genome}" \
        -O "${per_sample_gvcf_path}${sample}.g.vcf.gz" \
        -ERC GVCF

    # index GVCF file (skip if GATK already created it)
    if [[ ! -f "${per_sample_gvcf_path}${sample}.g.vcf.gz.tbi" ]]; then
        echo "Indexing GVCF for ${sample}..."
        tabix -p vcf "${per_sample_gvcf_path}${sample}.g.vcf.gz"
    fi

done

echo "Variant Calling Complete"
