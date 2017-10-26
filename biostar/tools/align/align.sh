set -ueo pipefail

# Get parameters.
INPUT={{data.path}}
GENOME={{genome.path}}
#ORGANISM={{organism.name}}

# Internal parameters.
RESULTS=results
DATA_DIR=$(dirname "$GENOME")
OUTPUT="mapped.bam"

# Build bwa index if not exist in project data folder.
if [ ! -f "$GENOME.bwt" ]; then
echo -e "Building bwa index.\n"
bwa index ${GENOME}
fi

# Align sequences
echo -e "Mapping reads to genome.\n"
mkdir -p ${RESULTS}
bwa mem -t 4 ${GENOME} ${INPUT} | samtools view -Sb |samtools sort >${RESULTS}/${OUTPUT}
samtools index ${RESULTS}/${OUTPUT}




