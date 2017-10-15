INPUT_DATA={{data.path}}

RESULT_DIR=results

mkdir -p $RESULT_DIR
fastqc -o $RESULT_DIR $INPUT_DATA
