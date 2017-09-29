BASE_DIR := test
RESDIR := test/res
WORKDIR := test/work
DATADIR := test/work/data
SCRIPTS := scripts
CHECK_INPUT := check_input.py
JOBID := job0


SAMPLEINFO := test/samples.txt
INPUTDATA := test/data.tar.gz
fname := $(notdir $(INPUTDATA))

untrimmed_out := $(WORKDIR)/multiqc_untrimmed/$(wildcard *.html)
all: setup_dir check_input $(untrimmed_out)


setup_dir:
	#
	# set up directories
	#
	mkdir -p $(RESDIR) $(WORKDIR)

extract_data:
	#
	# extract tar archive and move data to workdir
	tar -C $(WORKDIR) -xzvf data.tar.gz
	mkdir -p $(DATADIR)
	mv $(WORKDIR)/*.* $(DATADIR)


check_input:
	#
	# checks to see if all samples specified in the metadata are present in the input data.
	# And extracts them into $(WORKDIR)/data folder.
	#
	python $(SCRIPTS)/$(CHECK_INPUT) -s $(SAMPLEINFO) -i $(INPUTDATA) -w $(WORKDIR)

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

.PHONY: all clean
