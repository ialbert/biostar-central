
BATCH_SIZE=5000

# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

export POSTGRES_HOST=/var/run/postgresql

# Activate the conda environemnt.
conda activate engine

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

# Set the site domain.
SITE_DOMAIN=${1}
export SITE_DOMAIN=${SITE_DOMAIN}

# Add 5000 posts to search index every 3 minutes
python manage.py index --index ${BATCH_SIZE} --report