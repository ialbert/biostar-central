INPUT_DATA = {{data.path}}
INPUT_SAMPLE_INFO = {{sampleinfo.path}}

TRIM_QUALITY = {{quality_threshold.value}}
TRIM_PRIMER = {{ trim_primer.value }}
THREADS = {{ threads.value}}

DATA_DIR = data
WORK_DIR = work
RESULT_DIR  = results
ERR_LOG = stderr.txt

FASTQC_DIR =$(ANALYSIS_DIR)/fastqc
PTRIM_DIR = $(ANALYSIS_DIR)/trim_primer
QTRIM_DIR = $(ANALYSIS_DIR)/trim_quality

# Stores the updated sample information.
SAMPLE_INFO = updated_sampleinfo.txt


all:

	# Make a results directory.
	mkdir -p $(RESULT_DIR)

	# Extract reads from the archive.
	tar -xzvf $(INPUT_DATA)

	# Create a new samplesheet with the filenames.
	python -m biostar.tools.data.samplesheet $(INPUT_SAMPLE_INFO) $(DATA_DIR) > $(SAMPLE_INFO)

	# Create fastqc reports for all samples.
	mkdir -p $(FASTQC_DIR)
	cat $(SAMPLE_INFO) | parallel --verbose --progress -j $(THREADS) --header : --colsep '\t' fastqc  --nogroup -o $(FASTQC_DIR) $(DATA_DIR)/{file1} $(DATA_DIR)/{file2}

	# Run multiqc on the fastqc report.
	multiqc -f -n input_multiqc -o $(FASTQC_DIR) $(FASTQC_DIR)
	cp $(FASTQC_DIR)/input_multiqc.html $(RESULT_DIR)/initial_multiqc.html
	

{%  if trim_primer.value %}
	echo "Trimming primers"
	mkdir -p $(PTRIM_DIR)/logs

	# Run bbduk with primer trimming.
	cat $(SAMPLE_INFO) | parallel --verbose --progress -j $(THREADS)  --header : --colsep '\t' bbduk.sh \
	in1=$(DATA_DIR)/{file1} in2=$(DATA_DIR)/{file2} \
	out1=${sample_name}_R1_trimp.fq.gz \
	out2=$(PTRIM_DIR)/{sample_name}_R2_trimp.fq.gz \
	literal={fwd_primer},{rev_primer} ktrim=l k=15 hdist=1 tpe tbo  \
	overwrite=t stats=$(PTRIM_DIR)/logs/{sample_name}_trimp_stats.txt

	# Move statistic to results.
	cat $(PTRIM_DIR)/logs/*stats.txt > $(RESULT_DIR)/trimp_stats.txt

	# Run fastqc and
	ls $(PTRIM_DIR)/*gz | parallel --verbose --progress -j 8 fastqc --nogroup {} -o $(PTRIM_DIR)
	multiqc -n primmer_trimmed_multiqc --force -o $(PTRIM_DIR) $(PTRIM_DIR)
	cp $(PTRIM_DIR)/primmer_trimmed_multiqc.html $(RESULT_DIR)

    OUTPUT2 > INPUT
{% endif %}

{%  if trim_quality.value  %}



{% endif %}

 INPUT

	#
	# quality report
	#


#
# =============== Quality trimming ===================
# Trim reads to user specified threshold value.
#

# This here makes no sense


ifeq ($(wildcard $(PTRIM_DIR)),)
INDIR=$(DATA_DIR)
R1={file1}
R2={file2}
else
INDIR=$(PTRIM_DIR)
R1={sample_name}_R1_trimp.fq.gz
R2={sample_name}_R2_trimp.fq.gz
endif

trim_quality: $(RESULT_DIR) updated_sample_info
	echo $(INDIR),$(R1)
	cat  $(SAMPLE_INFO) |parallel --verbose --progress -j $(THREADS) --header : --colsep '\t' bbduk.sh \
	in1=$(INDIR)/$(R1) in2=$(INDIR)/$(R2) \
	qtrim=rl trimq=$(TRIM_QUALITY)  minlength=35 overwrite=true \
	out1=$(QTRIM_DIR)/{sample_name}_R1_trim.fq.gz \
	out2=$(QTRIM_DIR)/{sample_name}_R2_trim.fq.gz stats=$(QTRIM_DIR)/logs/{sample_name}_stats.txt
	cat $(QTRIM_DIR)/logs/*stats.txt > $(RESULT_DIR)/trimq_stats.txt
	#
	#
	# quality report
	#
	ls $(QTRIM_DIR)/*gz | parallel --verbose --progress -j 8 fastqc --nogroup {} -o $(QTRIM_DIR)
	multiqc -n quality_trimmed_multiqc --force -o $(QTRIM_DIR) $(QTRIM_DIR)
	cp $(QTRIM_DIR)/quality_trimmed_multiqc.html $(RESULT_DIR)

