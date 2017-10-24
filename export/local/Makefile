
testdata:
	mkdir -p lamar_data
	curl http://iris.bx.psu.edu/projects/metabarcode-data/data.tar.gz >lamar_data/data.tar.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/sampleinfo.txt >lamar_data/sampleinfo.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/1-SarriPal_S9_L001_R1_001.fastq.gz >lamar_data/1-SarriPal_S9_L001_R1_001.fastq.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/lamar_fishes_16S.fa >lamar_data/lamar_fishes_16S.fa
	curl http://iris.bx.psu.edu/projects/metabarcode-data/fish_accession_map.txt >lamar_data/fish_accession_map.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/accession_list.txt >lamar_data/accession_list.txt

centrifuge:
	conda install -y centrifuge

biopython:
	conda install -y biopython

taxonkit:
	conda install -y -c bioconda taxonkit

blast:
	# Install blast.
	conda install -y blast

nt_db:
	# Install nt database to export/local/nt
	# Ignores make error.
	mkdir -p nt
	(cd nt && update_blastdb.pl --decompress nt) || :

taxonomy:
	# Install taxonomy database to export/local/taxdb
	# Ignores make error.
	mkdir -p taxdb
	(cd taxdb && update_blastdb.pl --decompress taxdb) || :

blastdbs: blast taxonomy nt_db
