#
# A recipe demonstrating the use of parameters.
#

#
# Access parameters.
#
echo "Instrument: {{instrument.value}}"
echo "Protocol: {{protocol.value}}"
echo "Cutoff: {{cutoff.value}}"
echo "Validate: {{validate.value}}"
echo "Read length: {{readlen.value}}"

#
# The files to be used.
#
echo "Reads: {{reference.value}}"
echo "SRR: {{accession.value}}"

#
# Make a nested directory structure
#
mkdir -p data/store/

#
# Create a few files.
#
cat {{reference.value}} | grep  ">"  > data/store/sequence-names.txt
cat {{accession.value}} > data/store/srr.txt

#
# Generate a table of contents.
#
find . -name '*' > contents.txt

#
# Print the contents to the screen
#
cat contents.txt
