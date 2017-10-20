INPUT_DATA={{data.path}}

INPUT_FILE=$(basename "$INPUT_DATA")
SAMPLE="${INPUT_FILE%%.*}"

RESULT_DIR=results

mkdir -p $RESULT_DIR
fastqc -o $RESULT_DIR $INPUT_DATA


mv $RESULT_DIR/${SAMPLE}_fastqc.html {{settings.index}}
