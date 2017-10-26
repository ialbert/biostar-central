set -ueo pipefail

# Get parameters.
INPUT={{data.path}}
GENOME={{genome.path}}
ORGANISM={{organism.name}}

# Internal parameters.
RESULTS=results
DATA_DIR=$(dirname "$GENOME")
INDEX=${DATA_DIR}/index
OUTPUT="mapped.bam"

# Testing.

echo -e "input data: $INPUT"
echo -e "input genome: $GENOME"
echo -e "DATA_DIR: $DATA_DIR"
echo -e "INDEX: $INDEX"


# Build bwa index if not exist in project data folder.
if [ ! -f "$DATA_DIR/index.bwt" ]; then
echo -e "Building bwa index.\n"
bwa index -p index $GENOME
fi

# Align sequences
echo -e "Mapping reads to genome.\n"
bwa mem -t 4 $INDEX $INPUT | samtools view -Sb |samtools sort >${RESULTS}/$OUTPUT
samtools index ${RESULTS}/$OUTPUT




