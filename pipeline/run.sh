# run.sh $WORKDIR $SPEC $TEMPLATE $COMMAND

#set -ueo

WORKDIR=$1
SPEC=$2
TEMPLATE=$3
COMMAND=$4

COMMAND=${COMMAND:-all}

# check if inputs are valid

if [ ! -e "$SPEC" ] || [ ! -e  "$TEMPLATE" ] ; then
    echo "Input file not found; exiting"
    exit
fi

if [ ! -d "$WORKDIR" ]; then
    echo "$WORKDIR not found; Making"
    mkdir -p $WORKDIR
fi


# create make

python make.py $SPEC $TEMPLATE > $WORKDIR/Makefile
cd $WORKDIR

#exit
# run make
make $COMMAND

exit
#create report
#python metabarcode_results.py $report >index.html
