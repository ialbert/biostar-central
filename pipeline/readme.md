##  biostar-engine pipelines

### How to run
	
	make run
	  		
This will create ./test and runs everything in there.

To run in a different location

	make WORKDIR=./path/to/location

Currently it can run only inside pipeline.

### Variables in Makefile

	WORKDIR=./test
	SPEC=templates/metabarcode_qc/metabarcode_spec.hjson
	TEMPLATE=templates/metabarcode_qc/metabarcode_makefile.html
	COMMAND=all

### Make rules
	
	make workdir  ( creates work directory)
	make get_data ( get data and samplesheet )
	make setup    ( creates a Makefile using spec file and template)
	make run      ( runs all teh above commands)

### Results

Results will be in $WORKDIR/Results.

