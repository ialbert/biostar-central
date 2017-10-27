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
bwa mem -t 4 ${GENOME} ${INPUT} | samtools view -b |samtools sort >${OUTPUT}
samtools index ${OUTPUT}

# Get total no. of reads in input
TOTAL=$(bioawk -c fastx 'END{print NR}' ${INPUT})
MAPPED_READS=$(samtools view -F 4 ${OUTPUT} | cut -f 1 | sort |uniq |wc -l)

# Get number of reads mapped to each chromosome.
IDX_STATS=$(samtools idxstats ${OUTPUT})

# Collect results.
mkdir -p $RESULTS
mv ${OUTPUT}* ${RESULTS}/
echo -e "Total\tMapped" >${RESULTS}/mapping_stats.txt
echo -e "$TOTAL\t$MAPPED_READS" >>${RESULTS}/mapping_stats.txt
echo -e "$IDX_STATS" >${RESULTS}/chrom_mapping.txt
