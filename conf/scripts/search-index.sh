
BATCH_SIZE=5000

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

# Activate the conda environemnt.
conda activate engine

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

# Add 5000 posts to search index every 3 minutes
python manage.py index --index ${BATCH_SIZE} --report