# Data dump files.
DUMP_FILE=export/database/db.json

# Backup file.
BACKUP_DUMP_FILE=export/database/db.backup.`date +'%Y-%m-%d-%H%M'`.json

# Default settings module.
DJANGO_SETTINGS_MODULE := biostar.server.settings

# Default app.
DJANGO_APP :=

# Database name
DATABASE_NAME := database.db

# Command used to load initial data
LOAD_COMMAND := project

# Search index name
INDEX_NAME := index

# Search index directory
INDEX_DIR := search

# Recipes database to copy
COPY_DATABASE := recipes.db

all: recipes serve

accounts:
	$(eval DJANGO_SETTINGS_MODULE := biostar.accounts.settings)
	$(eval DJANGO_APP := biostar.accounts)

	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}

emailer:
	$(eval DJANGO_SETTINGS_MODULE := biostar.emailer.settings)
	$(eval DJANGO_APP := biostar.emailer)

	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}


pg:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}

recipes:
	$(eval DJANGO_SETTINGS_MODULE := biostar.recipes.settings)
	$(eval DJANGO_APP := biostar.recipes)
	$(eval LOAD_COMMAND := project)
	$(eval UWSGI_INI := site/test/recipes_uwsgi.ini)
	$(eval ANSIBLE_HOST := hosts/www.bioinformatics.recipes)
	$(eval ANSIBLE_ROOT := conf/ansible)
	$(eval SUPERVISOR_NAME := recipes)
	$(eval ENGINE_DIR := /export/www/biostar-central)

    # Set the settings variables.
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	@echo DATABASE_NAME=${DATABASE_NAME}


bioconductor:
	$(eval DJANGO_SETTINGS_MODULE := themes.bioconductor.settings)
	$(eval DJANGO_APP := biostar.forum)
	$(eval LOAD_COMMAND := populate)
	$(eval UWSGI_INI := themes/bioconductor/conf/uwsgi.ini)
	$(eval ANSIBLE_HOST := supportupgrade.bioconductor.org)
	$(eval ANSIBLE_ROOT := themes/bioconductor/conf/ansible)
	$(eval SUPERVISOR_NAME := forum)
	$(eval ENGINE_DIR := /export/www/biostar-central)

	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	@echo DATABASE_NAME=${DATABASE_NAME}


forum:
	$(eval DJANGO_SETTINGS_MODULE := biostar.forum.settings)
	$(eval DJANGO_APP := biostar.forum)
	$(eval LOAD_COMMAND := populate)
	$(eval UWSGI_INI := site/test/forum_uwsgi.ini)
	$(eval ANSIBLE_HOST := hosts/test.biostars.org)
	$(eval ANSIBLE_ROOT := conf/ansible)
	$(eval SUPERVISOR_NAME := engine)
	$(eval ENGINE_DIR := /export/www/biostar-engine)

	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	@echo DATABASE_NAME=${DATABASE_NAME}

serve:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py runserver --settings ${DJANGO_SETTINGS_MODULE}

init:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py collectstatic --noinput -v 0  --settings ${DJANGO_SETTINGS_MODULE}
	python manage.py migrate -v 0  --settings ${DJANGO_SETTINGS_MODULE}

load:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py loaddata --ignorenonexistent --settings ${DJANGO_SETTINGS_MODULE} $(DUMP_FILE)

delete:
	# Delete the database, logs and CACHE files.
	# Keep media and spooler.
	rm -rf export/logs/*.log
	rm -f export/db/${DATABASE_NAME}
	rm -rf export/static/CACHE
	rm -rf *.egg
	rm -rf *.egg-info

# Resets the site without removing jobs.
reset: delete init
    # Initializes the test project.

copy: reset
	@echo COPY_DATABASE=${COPY_DATABASE}
	python manage.py copy --db ${COPY_DATABASE} --settings ${DJANGO_SETTINGS_MODULE}

test:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}

	coverage run manage.py test ${DJANGO_APP} --settings biostar.server.test_settings -v 2 --failfast
	coverage html --skip-covered

	# Remove files associated with tests
	rm -rf export/tested


test_all:test

index:
	@echo INDEX_NAME=${INDEX_NAME}
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py index --settings ${DJANGO_SETTINGS_MODULE} --index 130000 --report

reindex:
	@echo INDEX_NAME=${INDEX_NAME}
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py index --remove --reset --index 3700000 --settings ${DJANGO_SETTINGS_MODULE}

demo: startup serve

startup:init
	python manage.py ${LOAD_COMMAND} --demo --settings ${DJANGO_SETTINGS_MODULE}

hard_reset: delete
	# Delete media and spooler.
	rm -rf export/spooler/*spool*
	rm -rf export/media/*

dump:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	python manage.py dumpdata --settings ${DJANGO_APP} --settings ${DJANGO_SETTINGS_MODULE} --exclude auth.permission --exclude contenttypes  > $(DUMP_FILE)
	@cp -f $(DUMP_FILE) $(BACKUP_DUMP_FILE)
	@ls -1 export/database/*.json

uwsgi: init
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo UWSGI_INI=${UWSGI_INI}
	uwsgi --ini ${UWSGI_INI}

transfer: pg
	python manage.py migrate --settings biostar.forum.settings
	python manage.py transfer -n 300 --settings biostar.transfer.settings

next:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py job --next --settings ${DJANGO_SETTINGS_MODULE}

config:
	(cd ${ANSIBLE_ROOT} && ansible-playbook -i ${ANSIBLE_HOST} config.yml --extra-vars -v)

install:
	(cd ${ANSIBLE_ROOT} && ansible-playbook -i ${ANSIBLE_HOST} install.yml --ask-become-pass --extra-vars -v)

deploy:
	(cd ${ANSIBLE_ROOT} && ansible-playbook -i ${ANSIBLE_HOST} server-deploy.yml --ask-become-pass --extra-vars "supervisor_program=${SUPERVISOR_NAME} restart=True engine_dir=${ENGINE_DIR}"  -v)
