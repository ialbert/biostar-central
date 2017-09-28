RESULTS_DIR=./res
WORKDIR=./work
DATADIR=./work/data
SCRIPTS=./scripts
CHECK_INPUT=check_input.py
JOBID=0

METADATA=sampleinfo.txt

INPUTDATA=data.tar.gz
INPUT_TYPE=gzipped_fastq



all: setup_dir check_input create_multiqc
        # execute all commands

setup_dir:
	#
	# set up directories
	#
	mkdir -p $(RESULTS_DIR) $(WORK_DIR)

check_input:
	#
	# checks to see if all samples specified in the metadata are present in the input data.
	# And extracts them into $(WORKDIR)/data folder.
	#
	python $(SCRIPTS)/$(CHECK_INPUT) -m $(METADATA) -i $(INPUTDATA) -w $(WORKDIR) --jobid $JOBID

QCDIR=$(WORKDIR)/multiQC
create_multiqc:
	#
	# run fastqc on all files
	#
	if [ -d "$(QCDIR)" ]; then rm -rf $(QCDIR); fi

	mkdir -p $(QCDIR)
	ls $(DATADIR)/*gz | parallel --verbose --progress -j 8 fastqc --nogroup {} -o $(QCDIR)
	#
	# run multiqc to create an aggregate report
	#
	multiqc -f -n multiqc $(QCDIR)
	cp multiqc.html $(RESULTS_DIR)
	mv multiqc* $(QCDIR)
	#
	# clean
	#
	rm -f $(QCDIR)/*fastqc.html $(QCDIR)/*fastqc.zip

