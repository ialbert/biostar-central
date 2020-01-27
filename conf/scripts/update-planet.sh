
UPDATE_COUNT=5


# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

# Activate the conda environemnt.
conda activate engine

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

# Update latest five entries for each planet blog.
python manage.py planet --update ${UPDATE_COUNT}
