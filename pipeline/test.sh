set -ue

WORKDIR=templates/metabarcode_qc
SPEC=templates/metabarcode_qc/metabarcode_spec.hjson
TEMPLATE=templates/metabarcode_qc/metabarcode_makefile.html


bash run.sh $WORKDIR $SPEC $TEMPLATE
