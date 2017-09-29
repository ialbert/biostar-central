WORKDIR=./work

JOBID := job0

fname := $(notdir $(INPUTDATA))
INPUTDATA := test/data.tar.gz

all : setup_dir extract_data fastqc multiqc

setup_dir:
	#
	# set up directories
	#
	mkdir -p $(WORKDIR)


extract_data:
	#
	# extract tar archive and move data to workdir
	tar -C $(WORKDIR) -xzvf $(INPUTDATA)
	#mkdir -p $(DATADIR)
	#mv $(WORKDIR)/*.* $(DATADIR)




MQCDIR := $(WORKDIR)/multiqc
QCDIR := $(WORKDIR)/fastqc
RESDIR := $(WORKDIR)/res

fastqc:
	#
	# run fastqc on all files
	#
	mkdir -p $(QCDIR)
	ls $(DATADIR)/*gz | parallel --verbose --progress -j 8 fastqc --nogroup {} -o $(QCDIR)

multiqc:
	#
	# run multiqc to create an aggregate report
	#
	if [ -d "$(MQCDIR)"   ]; then rm -rf $(MQCDIR); fi

	mkdir -p $(MQCDIR)

	multiqc -f -n multiqc $(QCDIR)
	cp multiqc.html $(RESDIR)
	mv multiqc* $(MQCDIR)
	#
	# clean
	#
	#rm -f $(QCDIR)/*fastqc.html
	rm -f $(QCDIR)/*fastqc.zip

 

