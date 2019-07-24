all: taxonomy input fishdb rnaseq

HOST=data.bioinformatics.recipes
USER=www
DIR=/export/sites/recipe_data/initial
DATA=recipes-initial-data.tar.gz

# Cleans initial data directories.
clean:
	rm -rf multiqc

# Packs up the current test data
data_pack:
	tar -cvzf ${DATA} multiqc

data_push:
	rsync -avz ${DATA} ${USER}@${HOST}:${DIR}

data_pull:
	rsync -avz ${USER}@${HOST}:${DIR}/${DATA} .
	tar zxvf ${DATA}

rnaseq: rnaseq/griffith/reads/HBR_1_R1.fq

rnaseq/griffith/reads/HBR_1_R1.fq:
	curl http://data.biostarhandbook.com/rnaseq/projects/griffith/griffith-data.tar.gz | tar xzv
	mkdir -p rnaseq/griffith/
	mv -f reads refs rnaseq/griffith/

# Taxonomy dependencies.
taxonomy: taxonomy/names.dmp taxonomy/nodes.dmp

taxonomy/names.dmp:
	centrifuge-download -o taxonomy taxonomy

taxonomy/nodes.dmp:
	centrifuge-download -o taxonomy taxonomy


# Test data.
input: input/fish_sequences.tar.gz

input/fish_sequences.tar.gz:
	mkdir -p input fish_species
	# scaffolds
	curl http://iris.bx.psu.edu/projects/testdata/assembly-data/scaffolds.fa >input/scaffolds.fasta
	# Get virus dataset
	curl http://iris.bx.psu.edu/projects/testdata/viral-data/col_bhv.fastq.gz > input/virus_sample1.fastq.gz
	curl http://iris.bx.psu.edu/projects/testdata/viral-data/D5.fastq.gz > input/virus_sample2.fastq.gz
	curl http://iris.bx.psu.edu/projects/testdata/viral-data/viral_genomes.fa >input/viral_genomes.fa
	# Get fish dataset
	curl http://iris.bx.psu.edu/projects/metabarcode-data/data.tar.gz  > input/fish_sequences.tar.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/sampleinfo.txt > input/lamar_sampleinfo.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/accession_list.txt > input/fish_accession_list.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/1-SarriPal_S9_L001_R1_001.fastq.gz > input/test_reads.fastq.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/1-SarriPal_S9_L001_R2_001.fastq.gz > input/test_reads_R2.fastq.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/fish_genomes.fa > fish_species/fish_species.fa
	cp input/test_reads.fastq.gz input/fish_metabarcode_sample1.fastq.gz
	cp input/test_reads.fastq.gz input/fish_metabarcode_sample2.fastq.gz
	cp input/test_reads.fastq.gz input/fish_metabarcode_reads.fastq.gz
	cp input/test_reads.fastq.gz input/test_reads_R1.fastq.gz


# Fish blast database.
fishdb: blastdb/fish_species/fish_species.nin

blastdb/fish_species/fish_species.nin: input
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



