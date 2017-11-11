set -ueo pipefail

# Locates all files in this directory
DATA_DIR=$(dirname {{sequence.path}})

# The output directory
OUTPUT=fastqc_results

# Create the output directory
mkdir -p $OUTPUT
mkdir -p $OUTPUT/zip

# Select all files in the directory that match fastq/fq
cat {{sequence.toc}} | egrep ".fastq|.fastq.gz|.fq|.fq.gz" | parallel -j 5 fastqc -o $OUTPUT {}

# Move the zip files out of sight
mv $OUTPUT/*.zip $OUTPUT/zip
