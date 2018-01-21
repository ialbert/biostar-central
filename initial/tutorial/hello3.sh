#
# Demonstrates the use of parameters.
#

# Print some of the parameters to the screen.
echo "Validate: {{validate.value}}"
echo "Speed: {{speed.value}}"
echo "Name: {{username.value}}"
echo "Max value: {{maxval.value}}"
echo "Direction: {{direction.value}}"
echo "Pace: {{pace.value}}"
echo "Reads: {{reads.value}}"
echo "Genes: {{genes.value}}"

#
# Make a nested directory structure
#
mkdir -p data/store/
mkdir -p data/temp/

#
# Create a files in directories.
#
echo "Direction={{direction.value}}" > data/store/direction.txt
echo "Pace={{pace.value}}" > data/store/pace.txt

#
# Generate a few test files.
#
ls -l > data/temp/ls.txt
find . -name '*' > data/temp/all.txt
