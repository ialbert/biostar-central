set -ueo pipefail

# Get parameters.
INPUT_DATA={{data.path}}
ACCESSION_LIST={{genome.path}}
ORGANISM_NAME={{organism.name}}

# Internal parameters.
REF_DIR=ref
RESULT_DIR=results

# Build bwa index.
bwa

