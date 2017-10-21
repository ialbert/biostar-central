INPUT_DATA={{data.path}}
INPUT_SAMPLE_INFO={{sampleinfo.path}}
TRIM_QUALITY={{quality_threshold.value}}
TRIM_PRIMER={{trim_primer.value}}
THREADS={{threads.value}}
KMER_LENGTH={{kmer_length.value}}
MIN_LENGTH={{read_length.value}}
RESULT_INDEX={{settings.index}}


# Set up files and folders.
DATA_DIR=data
RESULT_DIR=results
LOG=stdout.txt
FASTQC_DIR=work/fastqc


# Stores the updated sample information.
SAMPLE_INFO=updated_sampleinfo.txt

# Make a results directory.
mkdir -p ${RESULT_DIR}

# Extract reads from the archive.
tar -xzvf $INPUT_DATA

# Create a new samplesheet with the filenames.
python -m biostar.tools.data.samplesheet $INPUT_SAMPLE_INFO $DATA_DIR > $SAMPLE_INFO

# Sample sheet requires a work directory.
mkdir -p work

# Stores FASTQ results.
mkdir -p $FASTQC_DIR

# Create fastqc reports for all samples.
cat $SAMPLE_INFO | parallel --verbose --progress --header : --colsep '\t' fastqc  --nogroup -o $FASTQC_DIR {file1} {file2}

# Run multiqc on the fastqc report.
multiqc -f -n initial_multiqc -o work --no-data-dir $FASTQC_DIR
cp work/initial_multiqc.html ${RESULT_DIR}/initial_multiqc.html

# Copy files to input.
cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {file1} {result1}
cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {file2} {result2}

# Branching starts here.

{% if trim_primer.value %}

    # Trim primer.

    cat $SAMPLE_INFO | parallel --verbose --progress --header : --colsep '\t' bbduk.sh \
    in1={result1} in2={result2}  out1={temp1} out2={temp2} \
    literal={fwd_primer},{rev_primer} ktrim=l k=$KMER_LENGTH hdist=1 tpe tbo \
    overwrite=t stats=work/{sample_name}_primer_stats.txt

    # copy output to input
    cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {temp1} {result1}
    cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {temp2} {result2}

    # Move statistic to results.
    cat work/*primer_stats.txt > ${RESULT_DIR}/primer_trim_stats.txt

    # Remove previous files from $FASTQC_DIR.
    rm -f $FASTQC_DIR/*.zip $FASTQC_DIR/*.html

    # Run fastqc on all trimmed samples.
    cat $SAMPLE_INFO | parallel --verbose --progress --header : --colsep '\t' fastqc --nogroup -o $FASTQC_DIR {temp1} {temp2}

    # Create multiqc report of trimmed samples.
    multiqc -n primer_trimmed -o work --no-data-dir --force $FASTQC_DIR
    cp work/primer_trimmed.html ${RESULT_DIR}/primer_trimmed.html

{% endif %}

{% if trim_quality.value %}

    # Trim quality.

    cat  $SAMPLE_INFO |parallel --verbose --progress --header : --colsep '\t' bbduk.sh \
    in1={result1} in2={result2}  out1={temp1} out2={temp2} \
    qtrim=rl trimq=$TRIM_QUALITY  minlength=$MIN_LENGTH overwrite=true \
    stats=work/{sample_name}_qual_stats.txt

    # Copy output to input.
    cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {temp1} {result1}
    cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {temp2} {result2}

    # Move statistic to results.
    cat work/*qual_stats.txt > ${RESULT_DIR}/quality_trim_stats.txt

    # Remove previous files from $FASTQC_DIR
    rm -f $FASTQC_DIR/*.zip $FASTQC_DIR/*.html

    # Run fastqc on all trimmed samples.
    cat $SAMPLE_INFO | parallel --verbose --progress --header : --colsep '\t' fastqc --nogroup -o $FASTQC_DIR {temp1} {temp2}

    # Create multiqc report of trimmed samples.
    multiqc -n quality_trimmed -o work --no-data-dir --force $FASTQC_DIR
    cp work/quality_trimmed.html ${RESULT_DIR}/quality_trimmed.html

{% endif %}

# Copy trimmed data to results.
cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {result1} ${RESULT_DIR}/{sample_name}_trim_R1.fq.gz
cat $SAMPLE_INFO | parallel --header : --colsep '\t' cp {result2} ${RESULT_DIR}/{sample_name}_trim_R2.fq.gz

# collect results in $RESULT_INDEX
# TO DO
# python -m  biostar.tools.data.results $RESULT_INDEX template_file
cp $RESULT_DIR/quality_trimmed.html $RESULT_INDEX





