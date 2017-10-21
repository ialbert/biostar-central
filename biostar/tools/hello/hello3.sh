# Help the user track progress.
echo "Analysis started."

# Note how the variable matches the JSON structure.
FNAME={{foo.path}}

# Run a tool on the file.
ls -l $FNAME > result.txt

# All done.
echo "Analysis finished."