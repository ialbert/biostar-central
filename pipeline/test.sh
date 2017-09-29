analysis_name="metabarcode_qc"
spec="metabarcode_spec.hjson"
template="metabarcode_makefile.html"
report="metabarcode_report.html"

outdir=templates/${analysis_name}

# create make
python make.py  $outdir/$spec $outdir/$template >$outdir/Makefile

# run make
cd $outdir
make all

#create report
python metabarcode_results.py $report >index.html
