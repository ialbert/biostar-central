
all: taxonomy lamar genomes fishdb

blastdbs: taxdb nt

clean:
	rm -rf lamar taxonomy blastdb/fish genomes

# Taxonomy related files.
taxonomy/names.dmp:
	centrifuge-download -o taxonomy taxonomy

taxonomy/nodes.dmp:
	centrifuge-download -o taxonomy taxonomy

taxonomy: taxonomy/names.dmp taxonomy/nodes.dmp

# Lamar test data.
lamar/sequences.tar.gz:
	mkdir -p lamar
	curl http://iris.bx.psu.edu/projects/metabarcode-data/data.tar.gz  > lamar/sequences.tar.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/sampleinfo.txt > lamar/sampleinfo.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/accession_list.txt > lamar/accession_list.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/1-SarriPal_S9_L001_R1_001.fastq.gz > lamar/1-SarriPal_S9_L001_R1_001.fastq.gz


lamar: lamar/sequences.tar.gz

blastdb/fish/fish.fa: fish_genome
	mkdir -p blastdb/fish
	#
	# Get fish accession map.
	curl http://iris.bx.psu.edu/projects/metabarcode-data/fish_accession_map.txt >blastdb/taxid_map.txt
	#
	# Get acc2taxid table.
	#curl https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz | gunzip -c >blastdb/acc2taxid.txt
	#
	# Get taxid map in correct format.
	#cat blastdb/acc2taxid.txt | cut -f 2,3 | grep -v "accession" >blastdb/taxid_map.txt
	#
	# Make blast database.
	makeblastdb -in genomes/fish/fish.fa -dbtype nucl -out blastdb/fish/fish -parse_seqids -taxid_map blastdb/taxid_map.txt


# Fish blast database.
fishdb: blastdb/fish/fish.fa


# NCBI nt database.
blastdb/nt/*.gz:
	# Install nt database to export/local/blastdb/nt
	# Ignores make error.
	mkdir -p blastdb/nt
	(cd blastdb/nt && update_blastdb.pl --decompress nt) || :

nt:blastdb/nt/*.gz


# NCBI taxonomy database.
blastdb/taxdb/taxdb.bti:
	# Install taxonomy database to export/local/blastdb/taxdb
	# Ignores make error.
	mkdir -p blastdb/taxdb
	(cd blastdb/taxdb && update_blastdb.pl --decompress taxdb) || :

taxdb:blastdb/taxdb/taxdb.bti


# Fish genome.
genomes/fish/fish.fa:
	mkdir -p genomes/fish
	curl http://iris.bx.psu.edu/projects/metabarcode-data/fish_genomes.fa >genomes/fish/fish.fa

fish_genome:genomes/fish/fish.fa

genomes: fish_genome
