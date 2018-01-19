# Demonstrates the use of parameters

echo "Validate: {{validate.value}}"
echo "Speed: {{speed.value}}"
echo "Name: {{username.value}}"
echo "Max value: {{maxval.value}}"
echo "Direction: {{direction.value}}"
echo "Pace: {{pace.value}}"
echo "Reads: {{reads.value}}"
echo "Genes: {{genes.value}}"

# Make a nested directory structure
mkdir -p data/store/
mkdir -p data/temp/

# Create a files in directories.
echo "direction={{direction.value}} on reads={{reads.value}}" > data/store/{{direction.value}}.txt
echo "pace content {{pace.value}} on {{genes.value}}" > data/store/{{pace.value}}.txt
ls -l > data/temp/ls.txt
find . -name '*' > data/temp/all.txt
