set -ueo pipefail

#Get parameters.
INPUT_DATA={{sequence.path}}
INPUT_SAMPLE_INFO={{sampleinfo.path}}
THREADS={{threads.value}}
TRIM_QUALITY={{quality_threshold.value}}
KMER_LENGTH={{kmer_length.value}}
MIN_LENGTH={{read_length.value}}

# Main results.
RESULT_VIEW={{runtime.work_dir}}/{{settings.index}}

# Set up files and folders.
WORK_DIR={{runtime.work_dir}}/work
DATA_DIR={{runtime.work_dir}}/data
RESULT_DIR={{runtime.work_dir}}/results
FASTQC_DIR=${WORK_DIR}/fastqc

LC_ALL=C.UTF-8
LANG=C.UTF-8

# Stores the updated sample information.
SAMPLE_INFO=updated_sampleinfo.txt

# Make a results directory.
mkdir -p $RESULT_DIR

# Extract reads from the archive.
tar -xzvf $INPUT_DATA

# Create a new samplesheet with the filenames.
python -m biostar.tools.data.samplesheet $INPUT_SAMPLE_INFO $DATA_DIR > $SAMPLE_INFO

# Sample sheet requires a work directory.
mkdir -p $WORK_DIR

# Stores FASTQ results.
mkdir -p $FASTQC_DIR

# Create fastqc reports for all samples.
cat $SAMPLE_INFO | parallel --header : --colsep '\t' fastqc  --nogroup -o $FASTQC_DIR {file1} {file2}  2>/dev/null

# Run multiqc on the fastqc report.
echo -e "\n Creating multiqc report.\n"
multiqc -f -n initial_multiqc -o ${WORK_DIR} --no-data-dir $FASTQC_DIR 2>/dev/null
cp ${WORK_DIR}/initial_multiqc.html ${RESULT_DIR}/initial_multiqc.html


# Trim primer.

cat $SAMPLE_INFO | parallel  --header : --colsep '\t' bbduk.sh \
in1={file1} in2={file2}  out1={result1} out2={result2} \
literal={fwd_primer},{rev_primer} ktrim=l k=$KMER_LENGTH hdist=1 tpe tbo \
overwrite=t stats=${WORK_DIR}/{sample_name}_primer_stats.txt

# Move statistic to results.
cat ${WORK_DIR}/*primer_stats.txt > ${RESULT_DIR}/primer_trim_stats.txt

# Remove previous files from $FASTQC_DIR.
rm -f ${FASTQC_DIR}/*.zip ${FASTQC_DIR}/*.html

# Run fastqc on all trimmed samples.
cat $SAMPLE_INFO | parallel --header : --colsep '\t' fastqc --nogroup -o $FASTQC_DIR {result1} {result2}  2>/dev/null

# Create multiqc report of trimmed samples.
echo -e "\n Creating multiqc report.\n"
multiqc -n primer_trimmed -o ${WORK_DIR} --no-data-dir --force $FASTQC_DIR 2>/dev/null
cp ${WORK_DIR}/primer_trimmed.html ${RESULT_DIR}/primer_trimmed.html


# Trim quality.

cat  $SAMPLE_INFO |parallel --header : --colsep '\t' bbduk.sh \
in1={result1} in2={result2}  out1={temp1} out2={temp2} \
qtrim=rl trimq=$TRIM_QUALITY  minlength=$MIN_LENGTH overwrite=true \
stats=${WORK_DIR}/{sample_name}_qual_stats.txt

# Copy output to input.
cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {temp1} {result1}
cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {temp2} {result2}

# Move statistic to results.
cat ${WORK_DIR}/*qual_stats.txt > ${RESULT_DIR}/quality_trim_stats.txt

# Remove previous files from $FASTQC_DIR
rm -f ${FASTQC_DIR}/*.zip ${FASTQC_DIR}/*.html

# Run fastqc on all trimmed samples.
cat $SAMPLE_INFO | parallel --header : --colsep '\t' fastqc --nogroup -o $FASTQC_DIR {result1} {result2} 2>/dev/null

# Create multiqc report of trimmed samples.
echo -e "\n Creating multiqc report.\n"
multiqc -n quality_trimmed -o work --no-data-dir --force $FASTQC_DIR  2>/dev/null
cp ${WORK_DIR}/quality_trimmed.html ${RESULT_DIR}/quality_trimmed.html


# Collecting results.

# Copy trimmed data to results.
cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {result1} ${RESULT_DIR}/{sample_name}_trim_R1.fq.gz
cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {result2} ${RESULT_DIR}/{sample_name}_trim_R2.fq.gz

# Get read counts before trimming.
cat $SAMPLE_INFO | parallel --header : --colsep '\t' python -m biostar.tools.qc.qc_stats \
{sample_name} {file1}  >>${WORK_DIR}/raw_counts.txt

# Get read counts after trimming.
cat $SAMPLE_INFO | parallel --header : --colsep '\t' python -m biostar.tools.qc.qc_stats \
{sample_name} ${RESULT_DIR}/{sample_name}_trim_R1.fq.gz >>${WORK_DIR}/trimmed_counts.txt

# Get counts table.
cat ${WORK_DIR}/raw_counts.txt | sort -k 1,1 >${WORK_DIR}/raw_counts.s.txt
cat ${WORK_DIR}/trimmed_counts.txt | sort -k 1,1 >${WORK_DIR}/trimmed_counts.s.txt

echo -e "sample\traw\ttrimmed\n" >${RESULT_DIR}/count_stats.txt
join ${WORK_DIR}/raw_counts.s.txt  ${WORK_DIR}/trimmed_counts.s.txt >>${RESULT_DIR}/count_stats.txt

# Clean up.

rm -f ${WORK_DIR}/raw_counts*.txt  ${WORK_DIR}/trimmed_counts*.txt
rm -f ${WORK_DIR}/*.gz

# Copy main results to $RESULT_VIEW.

cp ${RESULT_DIR}/quality_trimmed.html $RESULT_VIEW
