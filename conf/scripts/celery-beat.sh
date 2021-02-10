# Load the conda commands.
source ~/miniconda3/etc/profile.d/conda.sh

export POSTGRES_HOST=/var/run/postgresql

# Activate the conda environemnt.
conda activate engine

# Set the configuration module.
export DJANGO_SETTINGS_MODULE=conf.run.site_settings

# Location of the log file
LOGFILE=/home/www/sites/biostar-central/live/logs/celery-beat.log

# The name of the application.
APP="biostar"

# The gunicorn instance to run.
CELERY="/home/www/.virtualenvs/biostar/bin/celery"

echo "starting celery beat with DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE"

$CELERY -A $APP beat -l info -f $LOGFILE

