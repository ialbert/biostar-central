set -ueo pipefail

# Get parameters.
INPUT={{data.path}}
GENOME={{genome.path}}

INDEX=${GENOME}
#ORGANISM={{organism.name}}
INDEX_DIR=$(dirname ${GENOME})/bwa

mkdir -p $INDEX_DIR
INDEX=${INDEX_DIR}/{{genome.uid}}

# Internal parameters.
WORK={{runtime.work_dir}}/work
mkdir -p $WORK

# Main results.
RESULT_VIEW={{runtime.work_dir}}/{{settings.index}}
OUTPUT=${WORK}/aligned.bam

# Build bwa index if not exist in project data folder.
if [ ! -f "$INDEX.bwt" ]; then
echo "Building bwa index."
bwa index -p ${INDEX} ${GENOME}
fi

# Align sequences
echo  "Mapping reads to the genome."
bwa mem -t 4 ${INDEX} ${INPUT}  | samtools view -b |samtools sort >${OUTPUT}
samtools index ${OUTPUT}

echo "Computing statistics."
echo "--------------------"
samtools flagstat ${OUTPUT}
echo "--------------------"

# Get statistics based on flags.
ALIGNMENTS=$(samtools view -c ${OUTPUT})
MAPPED=$(samtools view -c -F 4 ${OUTPUT})
MAPPED_FWD=$(samtools view -c -F 20 ${OUTPUT})
MAPPED_REV=$(samtools view -c -f 16 ${OUTPUT})
SECONDARY=$(samtools view -c -f 256 ${OUTPUT})
CHIMERIC=$(samtools view -c -f 2048 ${OUTPUT})
echo -e "Category\tCounts" > ${WORK}/alignment_stats.txt
echo -e "Mapped\t$MAPPED\nMapped_fwd\t$MAPPED_FWD\nMapped_rev\t$MAPPED_REV\n" >>${WORK}/alignment_stats.txt
echo -e "Secondary\t$SECONDARY\nChimeric\t$CHIMERIC\n" >>${WORK}/alignment_stats.txt


# Get number of reads mapped to each chromosome.
IDX_STATS=$(samtools idxstats ${OUTPUT})
echo -e "Chrom\tLength\tMapped\tUnmapped" >${WORK}/chrom_mapping.txt
echo -e "$IDX_STATS" >>${WORK}/chrom_mapping.txt


# Get total no. of reads in input
TOTAL=$(bioawk -c fastx 'END{print NR}' ${INPUT})
MAPPED_READS=$(samtools view -F 4 ${OUTPUT} | cut -f 1 | sort |uniq |wc -l)
UNMAPPED_READS=$(($TOTAL-$MAPPED_READS))


# Collect WORK.
echo -e "Total\tMapped\tUnmapped" >${WORK}/mapping_stats.txt
echo -e "$TOTAL\t$MAPPED_READS\t$UNMAPPED_READS" >>${WORK}/mapping_stats.txt

echo "Generating plots."

# Plot WORK.
python -m biostar.tools.align.plotter ${WORK}/chrom_mapping.txt ${WORK}/mapping_stats.txt >${RESULT_VIEW}

