# Frequently Asked Questions

### Contact

Contact email: [admin@biostars.org](mailto:admin@biostars.org)

### Interesting questions

* [What Are The Most Common Stupid Mistakes In Bioinformatics?](https://www.biostars.org/p/7126/)
* [Bioinformatics cartoon](https://www.biostars.org/p/16049/)

### Common questions

* How to convert Gene Symbols into Entrez ID and vice-varsa

Using [EntrezDirect](http://bit.ly/entrez-direct) 

```
# Remove `head -10` to get full list

$ esearch -db gene -query "human [orgn]" | efetch -format tabular | head -10
tax_id	Org_name	GeneID	CurrentID	Status	Symbol	Aliases	description	other_designations	map_location	chromosome	genomic_nucleotide_accession.version	start_position_on_the_genomic_accession	end_position_on_the_genomic_accession	orientation	exon_count	OMIM
9606	Homo sapiens	120893160	0	live	LOC120893160		Sharpr-MPRA regulatory region 11691		1q	1
9606	Homo sapiens	120893158	0	live	LOC120893158		Sharpr-MPRA regulatory region 11936		1q	1
9606	Homo sapiens	120893156	0	live	LOC120893156		Sharpr-MPRA regulatory region 5363		1p	1
9606	Homo sapiens	120893154	0	live	LOC120893154		Sharpr-MPRA regulatory region 8666		1p	1
9606	Homo sapiens	120893152	0	live	LOC120893152		Sharpr-MPRA regulatory region 9477		1p	1
9606	Homo sapiens	120893150	0	live	LOC120893150		Sharpr-MPRA regulatory region 3357		1p	1
9606	Homo sapiens	120893148	0	live	LOC120893148		Sharpr-MPRA regulatory region 2368		1p	1
9606	Homo sapiens	120893146	0	live	LOC120893146		Sharpr-MPRA regulatory region 5145		1p	1
9606	Homo sapiens	120893144	0	live	LOC120893144		Sharpr-MPRA regulatory region 4548		1p	1

```

### Support for Biostar

Biostar has been developed as an open source software with the **MIT licence** thanks to awards from the following
organizations:

* [US Fish and Wildlife Service Cooperative Agreement](https://www.fws.gov/grants/atc.html): award `F16AC01007`
* [National Institutes of Health (NIH)](http://www.nih.gov/), grant `NIH 5R25HG006243-02`
* [The Pennsylvania State University](http://www.psu.edu/)
