ASSEMBLY={{assembly.path}}
SPECIES={{species.value}}
NPROC={{processors.value}}

# Run augustus gene prediction.
echo "Running augustus and predcting genes."
echo "-------------------------------------"
cat $ASSEMBLY | parallel --j $NPROC --blocksize 5M --recstart '>' --pipe "cat {} > tmp/{%} && augustus --species=$SPECIES  tmp/{%}" >genes.gff
