
all: taxonomy lamar

# Taxonomy related files.
taxonomy/names.dmp:
	centrifuge-download -o taxonomy taxonomy

taxonomy/nodes.dmp:
	centrifuge-download -o taxonomy taxonomy

taxonomy: taxonomy/names.dmp taxonomy/nodes.dmp


# Lamar fish, test data.
lamar/data.tar.gz:
	mkdir -p lamar
	curl http://iris.bx.psu.edu/projects/metabarcode-data/data.tar.gz  > lamar/data.tar.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/sampleinfo.txt > lamar/sampleinfo.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/1-SarriPal_S9_L001_R1_001.fastq.gz > lamar/1-SarriPal_S9_L001_R1_001.fastq.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/lamar_fishes_16S.fa > lamar/lamar_fishes_16S.fa
	curl http://iris.bx.psu.edu/projects/metabarcode-data/fish_accession_map.txt > lamar/fish_accession_map.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/accession_list.txt > lamar/accession_list.txt

lamar: lamar/data.tar.gz


clean:
	rm -rf lamar taxonomy

blastdbs: blast tax_db nt_db

fishdb: blast fish_db blast_env

classify: centrifuge taxonomy

classify_tools: centrifuge taxonomy blast fish_db blast_env


nt_db:
	# Install nt database to export/local/blastdb/nt
	# Ignores make error.
	mkdir -p blastdb/nt
	mkdir -p blastdb/nt :

tax_db:
	# Install taxonomy database to export/local/blastdb/taxdb
	# Ignores make error.
	mkdir -p blastdb/taxdb
	(cd blastdb/taxdb && update_blastdb.pl --decompress taxdb) || :

fish_db:
	mkdir -p blastdb/fish
	#
	# Get fish sequences.
	curl http://iris.bx.psu.edu/projects/metabarcode-data/fish_genomes.fa >blastdb/fish/fish.fa
	curl http://iris.bx.psu.edu/projects/metabarcode-data/fish_accession_map.txt >blastdb/taxid_map.txt
	#
	# Get acc2taxid table.
	#curl https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz | gunzip -c >blastdb/acc2taxid.txt
	#
	# Get taxid map in correct format.
	#cat blastdb/acc2taxid.txt | cut -f 2,3 | grep -v "accession" >blastdb/taxid_map.txt
	#
	# make blastdb.
	makeblastdb -in blastdb/fish/fish.fa -dbtype nucl -out blastdb/fish/fish -parse_seqids -taxid_map blastdb/taxid_map.txt

blast_env:
	# set BLASTDB
	export BLASTDB=`pwd`/blastdb/fish









