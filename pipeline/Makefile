# set the working directory where the makefile will be created and analysis will be done.

WORKDIR=./test

# analysis-specific files.

SPEC=templates/metabarcode_qc/metabarcode_spec.hjson
TEMPLATE=templates/metabarcode_qc/metabarcode_makefile.html
COMMAND=all

# check if spec file exists.

ifeq ($(wildcard $(SPEC)),)
$(error $(SPEC) not found.)
endif

# check if template file exists.

ifeq ($(wildcard $(TEMPLATE)),)
$(error $(TEMPLATE) not found.)
endif


workdir:
	#
	# creates work directory.
	#
	mkdir -p $(WORKDIR)

get_data: workdir
	#
	# Download data and sampleinfo to workdir.
	#
	wget -P $(WORKDIR) http://iris.bx.psu.edu/projects/metabarcode-data/data.tar.gz
	wget -P $(WORKDIR) http://iris.bx.psu.edu/projects/metabarcode-data/sampleinfo.txt

setup:
	#
	# creates the makefile and runs the pipeline.
	#
	python make.py $(SPEC) $(TEMPLATE) > $(WORKDIR)/Makefile
	#
	# runs make
	#
	cd $(WORKDIR); make $(COMMAND)

run: workdir get_data setup
	@echo runs all commands to setup the pipeline and then runs it.
