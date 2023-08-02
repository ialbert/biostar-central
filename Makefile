
# Set the default Django settings module.
DJANGO_SETTINGS_MODULE ?= biostar.forum.settings

# Print usage information
usage:
	@echo "#"
	@echo "# Usage:"
	@echo "#"
	@echo "# 	make init        # initialize the database"
	@echo "# 	make serve       # run the development server"
	@echo "#"

# Initialize the database.
init:
	python manage.py collectstatic --noinput -v 0  --settings ${DJANGO_SETTINGS_MODULE}
	python manage.py migrate -v 0  --settings ${DJANGO_SETTINGS_MODULE}

# Development server
serve:
	python manage.py runserver --settings ${DJANGO_SETTINGS_MODULE}

# Run tests.
test:
	# Override the settings module.
	$(eval DJANGO_SETTINGS_MODULE=biostar.server.test_settings)
	coverage run manage.py test ${TEST} --settings ${DJANGO_SETTINGS_MODULE} -v 2 --failfast
	coverage html --skip-covered --omit="conf/*,biostar/celery.py,biostar/celeryconf.py"

	# Remove files associated with tests
	rm -rf export/tested

# Targets that do not correspond to a file.
.PHONY:  usage init serve test

