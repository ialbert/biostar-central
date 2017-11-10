set -ueo pipefail

# Get parameters.

GENOME={{genome.path}}
BAM={{bam.path}}

# Internal parameters.
GENOME_DIR=$(dirname ${GENOME})
GENOME_INDEX=${GENOME_DIR}/{{genome.uid}}.fai

BAM_DIR=$(dirname ${BAM})
BAM_INDEX=${BAM_DIR}/{{bam.uid}}.bai

# Result.
RESULT_VIEW={{runtime.work_dir}}/{{settings.index}}

# create genome if not exist.
if [ ! -f "$GENOME_INDEX" ]; then
echo "Building genome."
samtools faidx $GENOME
fi

# generate bam index if not exist.
if [ ! -f "$BAM_INDEX" ]; then
echo "Building bam index."
samtools index $BAM
fi

# create igv xml
python -m biostar.tools.igv.generate_xml $GENOME $BAM >$RESULT_VIEW
