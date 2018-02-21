# Stop script on any error.
set -ueo pipefail

# The directory that contains the data.
DATA_DIR=$(dirname {{sequence.value}})

# The output directory
OUTPUT=fastqc_results

# Create the output directory
mkdir -p $OUTPUT

# Select all files in the directory that match a known pattern.
cat {{sequence.toc}} | egrep ".fastq|.fastq.gz|.fq|.fq.gz|.bam" | parallel -j 5 fastqc -o $OUTPUT {}

# Move the zip files out of sight
mkdir -p $OUTPUT/zip
mv $OUTPUT/*.zip $OUTPUT/zip

# Generate a multiqc report.
# You may need to set the encoding in bash as
# export LC_ALL=C.UTF-8
# export LANG=C.UTF-8
multiqc fastqc_results -o fastqc_results
