
# Conda environment should be set

# Main codebase directory
SRC = ~/app/biostar-central

# The django settings module
SETTINGS = biostar.forum.settings

# Print usage information
usage:
	@echo "#"
	@echo "# Use the source, Luke!"
	@echo "#"

# Initialize the database.
init:
	cd ${SRC} && python manage.py migrate -v 0  --settings ${SETTINGS}
	python manage.py collectstatic --noinput -v 0  --settings ${SETTINGS}

# Development server
serve:
	python manage.py runserver --settings ${SETTINGS}

# Pull latest changes from github.
pull:
	cd ${SRC} && git pull

# Install requirements.
install:
	cd ${SRC} && pip install -r conf/requirements.txt

# Run tests.
test:
	# Override the settings module.
	$(eval DJANGO_SETTINGS_MODULE=biostar.server.test_settings)
	coverage run manage.py test ${TEST} --settings ${DJANGO_SETTINGS_MODULE} -v 2 --failfast
	coverage html --skip-covered --omit="conf/*,biostar/celery.py,biostar/celeryconf.py"

	# Remove files associated with tests
	rm -rf export/tested

# Various make related variables
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -c -u
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules


# Targets that do not correspond to a file.
.PHONY:  usage init serve migrate pull install test pg_data pg_drop pg_init pg_role pg_import pg_reset

