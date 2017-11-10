set -ueo pipefail

# The output directory
OUTPUT=fastqc
DATA_DIR=$(dirname {{sequence.path}})

cd {{runtime.work_dir}}

# Create the output directory
mkdir -p $OUTPUT


# Select all files in the data collection
find $DATA_DIR -name '*.fastq*' | parallel -j 5 fastqc -o $OUTPUT {}

# Select all files in the data collection
find $DATA_DIR -name '*fastq.qz*' | parallel -j 5 fastqc -o $OUTPUT {}
