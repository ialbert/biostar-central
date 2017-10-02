WORKDIR=./work
JOBID := job0

INPUTDATA := test/data.tar.gz
SAMPLEINFO := test/sampleinfo.txt
fname := $(notdir $(INPUTDATA))

PTRIM := True

all : check_input fastqc multiqc trim_primer trim_quality merge report

DATADIR:=$(WORKDIR)/data
RESDIR := $(WORKDIR)/res
check_input:
	python ../../check_input.py -i $(INPUTDATA) -s $(SAMPLEINFO) -w $(WORKDIR)
	mkdir -p $(RESDIR)


fastqc:
	#
	# run fastqc on all files
	#
	ls $(DATADIR)/*gz | parallel --verbose --progress -j 8 fastqc --nogroup {} -o $(DATADIR)

multiqc:
	#
	# run multiqc to create an aggregate report
	#
	multiqc -f -n multiqc -o $(DATADIR) $(DATADIR)
	cp $(DATADIR)/multiqc.html $(RESDIR)
	#
	# clean
	#
	#rm -f $(DATADIR)/*fastqc.zip



PTRIM_DIR=$(WORKDIR)/trim/trimp
trim_primer:
	#
	mkdir -p $(PTRIM_DIR)/logs
	#
	# trim left primer
	#
	cat $(SAMPLEINFO) |parallel --verbose --progress  --header : --colsep '\t' bbduk.sh \
	in1=$(DATADIR)/{file1} in2=$(DATADIR)/{file2} \
	out1=$(PTRIM_DIR)/{sample_name}_trimfwd1.fq.gz \
	out2=$(PTRIM_DIR)/{sample_name}_trimfwd2.fq.gz \
	literal={barcode}{fwd_primer} ktrim=l k=23 hdist=1 tpe tbo \
	stats=$(PTRIM_DIR)/{sample_name}_stats.txt overwrite=t skipr2=t

	#
	# trim right primer
	#
	cat $(SAMPLEINFO) | parallel --verbose --progress --header : --colsep '\t' bbduk.sh  \
	in1=$(PTRIM_DIR)/{sample_name}_trimfwd1.fq.gz \
	in2=$(PTRIM_DIR)/{sample_name}_trimfwd2.fq.gz \
	out1=$(PTRIM_DIR)/{sample_name}_R1_trim.fq.gz \
	out2=$(PTRIM_DIR)/{sample_name}_R2_trim.fq.gz  \
	literal={barcode}{rev_primer} ktrim=l k=28 hdist=1 tpe tbo \
	stats=$(PTRIM_DIR)/{sample_name}_stats.txt overwrite=t skipr1=t
	#
	# clean
	#
	mv $(PTRIM_DIR)/*stats.txt $(PTRIM_DIR)/logs
	rm -f $(PTRIM_DIR)/*trimfwd*



QTRIM_DIR=$(WORKDIR)/trim/trimq
trim_quality:
	mkdir -p $(QTRIM_DIR)/logs
ifeq ($(PTRIM),True)
	cat $(SAMPLEINFO) |parallel --verbose --progress  --header : --colsep '\t' bbduk.sh \
	in1=$(PTRIM_DIR)/{sample_name}_R1_trim.fq.gz in2=$(PTRIM_DIR)/{sample_name}_R2_trim.fq.gz \
	qtrim=rl trimq=30 overwrite=true minlength=35 \
	out1=$(QTRIM_DIR)/{sample_name}_R1_trimq.fq.gz \
	out2=$(QTRIM_DIR)/{sample_name}_R2_trimq.fq.gz
else
	cat $(SAMPLEINFO) |parallel --verbose --progress  --header : --colsep '\t' bbduk.sh \
	in1=$(DATADIR)/{file1} in2=$(DATADIR)/{file2} \
	qtrim=rl trimq=30 overwrite=true minlength=35 \
	out1=$(QTRIM_DIR)/{sample_name}_R1_trimq.fq.gz
	out2=$(QTRIM_DIR)/{sample_name}_R2_trimq.fq.gz
endif
	#
	# quality report
	#
	ls $(QTRIM_DIR)/*gz | parallel --verbose --progress -j 8 fastqc --nogroup {} -o $(QTRIM_DIR)
	multiqc -n qtrim_multiqc -o $(QTRIM_DIR) $(QTRIM_DIR)
	cp $(QTRIM_DIR)/qtrim_multiqc.html $(RESDIR)



MERGE_DIR=$(WORKDIR)/merge
merge:
	mkdir -p $(MERGE_DIR)/logs
	if [ -d $(QTRIM_DIR) ]; then \
		cat $(SAMPLEINFO) |parallel --verbose --progress  --header : --colsep '\t' bbmerge.sh \
		in1=$(QTRIM_DIR)/{sample_name}_R1_trimq.fq.gz in2=$(QTRIM_DIR)/{sample_name}_R2_trimq.fq.gz \
		out=$(MERGE_DIR)/{sample_name}_merged.fq.gz \
		mismatches=3 minoverlap=12 ; \
	elif [ -d $(PTRIM_DIR) ] ; then  \
		cat $(SAMPLEINFO) |parallel --verbose --progress  --header : --colsep '\t' bbmerge.sh \
		in1=$(PTRIM_DIR)/{sample_name}_R1_trim.fq.gz in2=$(PTRIM_DIR)/{sample_name}_R2_trim.fq.gz \
		out=$(MERGE_DIR)/{sample_name}_merged.fq.gz \
		mismatches=3 minoverlap=12 ; \
	else \
		cat $(SAMPLEINFO) |parallel --verbose --progress  --header : --colsep '\t' bbmerge.sh \
		in1=$(DATADIR)/{sample_name}_R1_trim.fq.gz in2=$(DATADIR)/{sample_name}_R2_trim.fq.gz \
		out=$(MERGE_DIR)/{sample_name}_merged.fq.gz \
		mismatches=3 minoverlap=12 ; \
	fi ;
	#
	# quality report
	#
	ls $(MERGE_DIR)/*gz | parallel --verbose --progress -j 8 fastqc --nogroup {} -o $(MERGE_DIR)
	multiqc -n merge_multiqc -o $(MERGE_DIR) $(MERGE_DIR)
	cp $(MERGE_DIR)/merge_multiqc.html $(RESDIR)


report:
	#
	# create html report of result.
	# inputs are result dir and template
	#
	python metabarcode_results.py $(RESDIR) metabarcode_report.html > index.html







