READ1={{read1.path}}
READ2={{read2.path}}
MAX_LEN={{maxlength.value}}
RESULT_VIEW={{runtime.work_dir}}/{{settings.index}}

# Internal parameters.
RESULT_DIR=results
OUT=merged.fq.gz
FILTERED=mergedWithoutNs.fq.gz

# Make $RESULT_DIR
mkdir -p $RESULT_DIR

# Merge read pairs.
echo "Merging reads."
echo "--------------------"
bbmerge.sh in1=$READ1 in2=$READ2 out=${RESULT_DIR}/$OUT maxlength=$MAX_LEN 2>&1

# FastQC report
fastqc --nogroup ${RESULT_DIR}/$OUT 2>/dev/null

# Move fastqc report to main result view.
mv ${RESULT_DIR}/*.html $RESULT_VIEW
rm -f ${RESULT_DIR}/*fastqc.zip

# Remove reads with Ns.
echo "--------------------------"
echo "Filtering reads with Ns."
echo "--------------------------"
bbduk.sh in=${RESULT_DIR}/$OUT out=${RESULT_DIR}/$FILTERED maxns=0 2>&1



