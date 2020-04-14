# Data dump files.
DUMP_FILE=export/database/db.json

# Backup file.
BACKUP_DUMP_FILE=export/database/db.backup.`date +'%Y-%m-%d-%H%M'`.json

# Default settings module.
DJANGO_SETTINGS_MODULE := themes.bioconductor.settings

# Default app.
DJANGO_APP := biostar.forum

# Database name
DATABASE_NAME := bioconductor.db

# Search index name
INDEX_NAME := index

# Search index directory
INDEX_DIR := search

$(eval UWSGI_INI := themes/bioconductor/conf/uwsgi.ini)
$(eval ANSIBLE_HOST := supportupgrade.bioconductor.org)
$(eval ANSIBLE_ROOT := themes/bioconductor/conf/ansible)
$(eval SUPERVISOR_NAME := forum)
$(eval ENGINE_DIR := /export/www/biostar-central)

@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
@echo DJANGO_APP=${DJANGO_APP}
@echo DATABASE_NAME=${DATABASE_NAME}

all: recipes serve

pg:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}

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
	python manage.py populate --demo --settings ${DJANGO_SETTINGS_MODULE}

dump:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	python manage.py dumpdata --settings ${DJANGO_SETTINGS_MODULE} --exclude auth.permission --exclude contenttypes  > $(DUMP_FILE)
	@cp -f $(DUMP_FILE) $(BACKUP_DUMP_FILE)
	@ls -1 export/database/*.json

uwsgi: init
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo UWSGI_INI=${UWSGI_INI}
	uwsgi --ini ${UWSGI_INI}

config:
	ansible-playbook -i ${HOST} ../conf/ansible/server-config.yml -v

install:
	ansible-playbook -i ${HOST} ../conf/ansible/server-install.yml -v

deploy:
	ansible-playbook -i ${HOST} ../conf/ansible/server-deploy.yml --ask-become-pass --extra-vars "restart=True start_biocond=True"  -v
