
INPUT_DATA={{data.path}}
TAXA_MAP={{taxamap.path}}
INDEX=${INPUT_DATA%.*}

# Download taxonomy
centrifuge-download -o taxonomy taxonomy

# Build centrifuge index.
centrifuge-build -p 4 --conversion-table $TAXA_MAP --taxonomy-tree taxonomy/nodes.dmp \
--name-table taxonomy/names.dmp $INPUT_DATA $INDEX