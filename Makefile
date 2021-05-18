# Database JSON dump files.
SAVE_FILE=export/backup/db.last.json

# Backup file.
BACKUP_FILE=export/backup/db.`date +'%Y-%m-%d-%H%M'`.json

# Default settings module.
DJANGO_SETTINGS_MODULE := biostar.server.settings

# Default app.
DJANGO_APP :=

# Command used to load initial data
LOAD_COMMAND := project

# Search index name
INDEX_NAME := index

# What test to run, defaults to all tests
TEST := biostar

# Search index directory
INDEX_DIR := search
# Some variables need to come from the enviroment.

.PHONY:  recipes accounts demo

# Recipes database to copy
export COPY_DATABASE := recipes.db

# Database name is accessed via an enviroment variable.
export DATABASE_NAME := database.db

all: recipes serve

accounts:
	@echo "*** Setting variables for accounts app."
	$(eval DJANGO_SETTINGS_MODULE := biostar.accounts.settings)
	$(eval CELERY_CONF := biostar/celeryconf.py)
	$(eval DJANGO_APP := biostar.accounts)
	$(eval TEST:=biostar.accounts)

emailer:
	@echo "*** Setting variables for emails app."
	$(eval DJANGO_SETTINGS_MODULE := biostar.emailer.settings)
	$(eval DJANGO_APP := biostar.emailer)
	$(eval TEST:=biostar.emailer)

recipes:
	@echo "*** Setting variables for recipe app."
	$(eval DJANGO_SETTINGS_MODULE := biostar.recipes.settings)
	$(eval DJANGO_APP := biostar.recipes)
	$(eval LOAD_COMMAND := project)
	$(eval UWSGI_INI := conf/site/site_uwsgi.ini)
	$(eval WSGI_FILE := biostar/recipes/wsgi.py)
	$(eval TASKS_MODULE := biostar.recipes.tasks)
	$(eval TARGET:=recipes)
	$(eval TEST:=biostar.recipes)


bioconductor:
	@echo "*** Setting variables for bioconductor app."
	$(eval DJANGO_SETTINGS_MODULE := themes.bioconductor.settings)
	$(eval DJANGO_APP := biostar.forum)
	$(eval LOAD_COMMAND := populate)
	$(eval UWSGI_INI := conf/site/site_uwsgi.ini)
	$(eval WSGI_FILE := themes/bioconductor/wsgi.py)
	$(eval TASKS_MODULE := biostar.forum.tasks)
	$(eval TARGET:=supportupgrade)
	$(eval TEST:=biostar.forum)

forum:
	@echo "*** Setting variables for forum app."
	$(eval DJANGO_SETTINGS_MODULE := biostar.forum.settings)
	$(eval DJANGO_APP := biostar.forum)
	$(eval LOAD_COMMAND := populate)
	$(eval UWSGI_INI := conf/site/site_uwsgi.ini)
	$(eval TASKS_MODULE := biostar.forum.tasks)
	$(eval WSGI_FILE := biostar/forum/wsgi.py)
	$(eval TEST:=biostar.forum)

biostar: forum
	$(eval TARGET := biostar)

test:
	$(eval TARGET := test)

echo:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	@echo DATABASE_NAME=${DATABASE_NAME}

serve: init
	python manage.py runserver --settings ${DJANGO_SETTINGS_MODULE}

init: echo
	python manage.py collectstatic --noinput -v 0  --settings ${DJANGO_SETTINGS_MODULE}
	python manage.py migrate -v 0  --settings ${DJANGO_SETTINGS_MODULE}

runtest:
	@echo DJANGO_SETTINGS_MODULE=biostar.server.test_settings
	@echo DJANGO_APP=${DJANGO_APP}
	$(eval DJANGO_SETTINGS_MODULE=biostar.server.test_settings)
	coverage run manage.py test ${TEST} --settings ${DJANGO_SETTINGS_MODULE} -v 2 --failfast
	coverage html --skip-covered --omit="conf/*,biostar/celery.py,biostar/celeryconf.py"

	# Remove files associated with tests
	rm -rf export/tested

test_all:runtest

index:
	@echo INDEX_NAME=${INDEX_NAME}
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py index --settings ${DJANGO_SETTINGS_MODULE} --index 130000 --report

reindex:
	@echo INDEX_NAME=${INDEX_NAME}
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py index --remove --reset --index 3700000 --settings ${DJANGO_SETTINGS_MODULE}

demo: startup serve

startup: init
	python manage.py ${LOAD_COMMAND} --demo --settings ${DJANGO_SETTINGS_MODULE}

copy: reset
	@echo COPY_DATABASE=${COPY_DATABASE}
	python manage.py copy --db ${COPY_DATABASE} --settings ${DJANGO_SETTINGS_MODULE}

reset: echo
	# Delete the database, logs and CACHE files.
	# Keep media and spooler.
	rm -rf export/logs/*.log
	rm -rf export/spammers/
	# Database is always found in export/db/
	rm -f export/db/${DATABASE_NAME}
	rm -rf export/static/CACHE
	rm -rf *.egg
	rm -rf *.egg-info

hard_reset: reset
	# Delete media and spooler.
	rm -rf export/spooler/*spool*
	rm -rf export/media/*

load:
	# Loads a data fixture.
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py loaddata --ignorenonexistent --settings ${DJANGO_SETTINGS_MODULE} $(SAVE_FILE)

save:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	python manage.py dumpdata ${DJANGO_APP} --settings ${DJANGO_SETTINGS_MODULE} --exclude auth.permission --exclude contenttypes  > $(SAVE_FILE)
	@cp -f $(SAVE_FILE) $(BACKUP_FILE)
	@ls -1 export/backup/*.json

uwsgi: init
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo WSGI_FILE=${WSGI_FILE}
	uwsgi --http :8000 --wsgi-file ${WSGI_FILE} --master --spooler export/spooler/ --pythonpath $(shell pwd) --import ${TASKS_MODULE}

transfer:
	python manage.py migrate --settings biostar.forum.settings
	python manage.py transfer -n 300 --settings biostar.transfer.settings

next:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py job --next --settings ${DJANGO_SETTINGS_MODULE}

config:
	cd conf/ansible && make config

install:
	cd conf/ansible && make install

redis:
	# Turn on redis
	redis-server

celery:
	# Run celery ( including beat )
	# Requires redis to be ran in a different shell.
	@echo CELERY_CONF=${CELERY_CONF}
	celery -A biostar worker -B -l INFO

remote_transfer:
	cd conf/ansible && make transfer

deploy:
	@echo TARGET=${TARGET}
	(cd conf/ansible && make deploy TARGET=${TARGET})

# Temporary command to deploy forum on test server - with directory name being ~/biostar-engine/
# Being refactored out once migrated.
forum_deploy:
	@echo "Deploying forum on test server."
	(cd conf/ansible && make forum_deploy TARGET=test)

clear_cache:
	echo 'flush_all' | nc localhost 11211