
INPUT_DATA={{data.path}}
TAXA_MAP={{taxamap.path}}
REFERENCE={{reference.path}}
INDEX=${REFERENCE%.*}
TAX_NAMES={{tax_names.path}}
TAX_NODES={{tax_nodes.path}}

DATA_DIR=data
RESULT_DIR=results

# Build centrifuge index.
#centrifuge-build -p 4 --conversion-table $TAXA_MAP --taxonomy-tree $TAX_NODES \
#--name-table $TAX_NAMES $REFERENCE $INDEX

# Run classification using centrifuge.
mkdir -p $RESULT_DIR
centrifuge -x $INDEX -U $INPUT_DATA --report-file ${RESULT_DIR}/report.txt -S ${RESULT_DIR}/classification.txt

# Register results to result/index.html
# TO DO



