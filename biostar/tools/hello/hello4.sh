# Help the user track progress.
echo "Analysis started."

# Run a simple unix tool on the file.
ls -l {{data.value}} > results.txt

# All done.
echo "Analysis finished."
