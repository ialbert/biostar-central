all: taxonomy lamar fishdb

clean:
	rm -rf blastdb lamar taxonomy input fish_species

# Taxonomy dependencies.
taxonomy: taxonomy/names.dmp taxonomy/nodes.dmp

taxonomy/names.dmp:
	centrifuge-download -o taxonomy taxonomy

taxonomy/nodes.dmp:
	centrifuge-download -o taxonomy taxonomy


# Lamar test data.
lamar: lamar/sequences.tar.gz

lamar/sequences.tar.gz:
	mkdir -p input fish_species
	curl http://iris.bx.psu.edu/projects/metabarcode-data/data.tar.gz  > input/lamar_sequences.tar.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/sampleinfo.txt > input/lamar_sampleinfo.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/accession_list.txt > input/accession_list.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/1-SarriPal_S9_L001_R1_001.fastq.gz > input/fastqc_test_reads.fastq.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/fish_genomes.fa > fish_species/fish_species.fa
	cp input/fastqc_test_reads.fastq.gz input/align_test_reads.fastq.gz
	cp input/fastqc_test_reads.fastq.gz input/classify_test_fish_metabarcode_reads.fastq.gz
	cp input/fastqc_test_reads.fastq.gz input/align_classify_test_metabarcode_reads.fastq.gz

# Fish blast database.
fishdb: blastdb/fish_species/fish_species.nin

blastdb/fish_species/fish_species.nin: lamar
	mkdir -p blastdb/fish_species
	#
	# Get fish accession map
	curl http://iris.bx.psu.edu/projects/metabarcode-data/fish_accession_map.txt > blastdb/taxid_map.txt
	#
	# Get acc2taxid table.
	#curl https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz | gunzip -c > blastdb/acc2taxid.txt
	#
	# Get taxid map in correct format.
	#cat blastdb/acc2taxid.txt | cut -f 2,3 | grep -v "accession" > blastdb/taxid_map.txt
	#
	# Make blast database.
	makeblastdb -in fish_species/fish_species.fa -dbtype nucl -out blastdb/fish_species/fish_species -parse_seqids -taxid_map blastdb/taxid_map.txt


# Blast  public databases.
blastdbs: taxdb nt

# NCBI nt database.
nt:nt/*.gz

nt/*.gz:
	# Install nt database to export/local/blastdb/nt
	# Ignores make error.
	mkdir -p nt
	(cd nt && update_blastdb.pl --decompress nt) || :


# NCBI taxonomy database.
taxdb:taxdb/taxdb.bti

taxdb/taxdb.bti:
	# Install taxonomy database to export/local/blastdb/taxdb
	# Ignores make error.
	mkdir -p taxdb
	(cd taxdb && update_blastdb.pl --decompress taxdb) || :



