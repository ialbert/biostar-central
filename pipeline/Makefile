# set the working directory where the makefile will be created and analysis will be done.

TEMP_DIR = ~/tmp

# analysis-specific files.

SPEC=templates/qc/qc_spec.hjson
TEMPLATE=templates/qc/qc_makefile.html
COMMAND=all

TEST_DATA=$(TEMP_DIR)/data.tar.gz
SAMPLE_DATA=$(TEMP_DIR)/sampleinfo.txt

run: setup
	@echo runs all commands to setup the pipeline and then runs it.

$(TEMP_DIR):
	mkdir -p $(TEMP_DIR)

$(TEST_DATA): $(TEMP_DIR)
	curl -o $(TEST_DATA) -O http://iris.bx.psu.edu/projects/metabarcode-data/data.tar.gz

$(SAMPLE_DATA): $(TEMP_DIR)
	 curl -o $(SAMPLE_DATA) -O http://iris.bx.psu.edu/projects/metabarcode-data/sampleinfo.txt

setup: $(TEST_DATA) $(SAMPLE_DATA)
	#
	# creates the makefile and runs the pipeline.
	#
	#python make.py $(SPEC) $(TEMPLATE) > $(WORKDIR)/Makefile
	#
	# runs make
	#
	#cd $(WORKDIR); make $(COMMAND)

