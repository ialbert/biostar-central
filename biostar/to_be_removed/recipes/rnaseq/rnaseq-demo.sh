#!/usr/bin/env bash

# Exit this script on any error.
set -euxo pipefail

# This is the directory that stored the reads
DATA={{reads.data_dir}}

# Set the reference genome.
REF={{reference.value}}

# The annotations.
GTF={{annotation.value}}

# Collect the output of more verbose commands here.
RUNLOG=runlog.txt

# Create output folder for BAM files.
mkdir -p bam

# The index determines what the data is aligned against.
mkdir -p index
IDX=index/reference

# Build the hisat2 index.
hisat2-build $REF $IDX 1> $RUNLOG 2> $RUNLOG

# Align the files.
for SAMPLE in HBR UHR;
do
    for REPLICATE in 1 2 3;
    do
        # Build the name of the files.
        R1=${DATA}/${SAMPLE}_${REPLICATE}_R1.fq
        R2=${DATA}/${SAMPLE}_${REPLICATE}_R2.fq

        # The name of the BAM files.
        BAM=bam/${SAMPLE}_${REPLICATE}.bam

        # Run the aligner.
        echo "*** Aligning: $BAM"
        hisat2 $IDX -1 $R1 -2 $R2 2>> $RUNLOG | samtools sort > $BAM 2>> $RUNLOG
        samtools index $BAM
    done
done

# Perform the feature counting
echo "*** Counting features with: $GTF"
featureCounts -a $GTF -g gene_name -o counts.txt  bam/HBR*.bam  bam/UHR*.bam 2>> $RUNLOG

echo "*** Generating simple counts."
cat counts.txt | cut -f 1,7-12 > simple_counts.txt

#
# We run the data through three methods: deseq1, deseq2 and edgeR
# For each we produce the enrichment table and a heatmap.
#

#
# Get the R scripts from the Handbook
#
# This directory will store the scripts
mkdir -p scripts

curl -s http://data.biostarhandbook.com/rnaseq/code/deseq1.r > scripts/deseq1.r
curl -s http://data.biostarhandbook.com/rnaseq/code/deseq2.r > scripts/deseq2.r
curl -s http://data.biostarhandbook.com/rnaseq/code/edger.r > scripts/edger.r
curl -s http://data.biostarhandbook.com/rnaseq/code/draw-heatmap.r > scripts/draw-heatmap.r

# Collect the results in this folder.
mkdir -p results

#
# Run each of the three methods
#

#
# Deseq1
#
echo "*** Running Deseq1"
cat simple_counts.txt | Rscript scripts/deseq1.r 3x3 > results/deseq1.txt  2>> $RUNLOG
cat norm-matrix-deseq1.txt | Rscript scripts/draw-heatmap.r > results/deseq1_heatmap.pdf

#
# DESeq2
#
echo "*** Running Deseq2"
cat simple_counts.txt | Rscript scripts/deseq2.r 3x3 > results/deseq2.txt  2>> $RUNLOG
cat norm-matrix-deseq2.txt | Rscript scripts/draw-heatmap.r > results/deseq2_heatmap.pdf

#
# EdgeR
#
echo "*** Running edgeR"
cat simple_counts.txt | Rscript scripts/edger.r 3x3 > results/edger.txt  2>> $RUNLOG
cat norm-matrix-edgeR.txt | Rscript scripts/draw-heatmap.r > results/edger_heatmap.pdf
