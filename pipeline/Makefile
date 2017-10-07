TEMP_DIR = ~/tmp/qc_test
SAMPLE_DATA = $(TEMP_DIR)/sampleinfo.txt
TEST_DATA = $(TEMP_DIR)/data.tar.gz

RUNDIR=~/tmp


all: run
	@echo runs all commands to setup the pipeline and then runs it.

$(TEMP_DIR):
	mkdir -p $(TEMP_DIR)

$(TEST_DATA): $(TEMP_DIR)
	curl  http://iris.bx.psu.edu/projects/metabarcode-data/data.tar.gz > $(TEST_DATA)
	#(cd $(TEMP_DIR) && tar -xzvmf data.tar.gz)

$(SAMPLE_DATA): $(TEMP_DIR)
	 curl http://iris.bx.psu.edu/projects/metabarcode-data/sampleinfo.txt > $(SAMPLE_DATA)

getdata: $(TEST_DATA) $(SAMPLE_DATA)

run: getdata
	python -m pipeline.qc.test > $(RUNDIR)/Makefile
	cd $(RUNDIR); make
