# Data dump files.
DUMP_FILE=../../export/database/db.json

# Backup file.
BACKUP_DUMP_FILE=../../export/database/db.backup.`date +'%Y-%m-%d-%H%M'`.json

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


all: serve


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


serve:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	@echo UWSGI_INI=${UWSGI_INI}
	@echo DATABASE_NAME=${DATABASE_NAME}
	@echo ANSIBLE_HOST=${ANSIBLE_HOST}
	@echo ANSIBLE_ROOT=${ANSIBLE_ROOT}
	python ../../manage.py runserver --settings ${DJANGO_SETTINGS_MODULE}


init:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python ../../manage.py collectstatic --noinput -v 0  --settings ${DJANGO_SETTINGS_MODULE}
	python ../../manage.py migrate -v 0  --settings ${DJANGO_SETTINGS_MODULE}

load:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python ../../manage.py loaddata --ignorenonexistent --settings ${DJANGO_SETTINGS_MODULE} $(DUMP_FILE)

delete:
	# Delete the database, logs and CACHE files.
	# Keep media and spooler.
	rm -rf ../../export/logs/*.log
	rm -f ../../export/db/${DATABASE_NAME}
	rm -rf ../../export/static/CACHE
	rm -rf ../../*.egg
	rm -rf ../../*.egg-info

# Resets the site without removing jobs.
reset: delete init
    # Initializes the test project.

copy: reset
	@echo COPY_DATABASE=${COPY_DATABASE}
	python manage.py copy --db ${COPY_DATABASE} --settings ${DJANGO_SETTINGS_MODULE}

test:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	#python manage.py test ${DJANGO_APP} --settings ${DJANGO_SETTINGS_MODULE} -v 2 --failfast
	coverage run ../../manage.py test ${DJANGO_APP} --settings ${DJANGO_SETTINGS_MODULE} -v 2 --failfast
	coverage html --skip-covered

test_all:
	python ../../manage.py test --settings biostar.test.test_settings -v 2 --failfast
	coverage run ../../manage.py test --settings biostar.test.test_settings -v 2 --failfast
	coverage html --skip-covered

index:
	@echo INDEX_NAME=${INDEX_NAME}
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python manage.py index --settings ${DJANGO_SETTINGS_MODULE} --index 130000 --report

reindex:
	@echo INDEX_NAME=${INDEX_NAME}
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	python ../../manage.py index --remove --reset --index 130000 --settings ${DJANGO_SETTINGS_MODULE}

forum_demo: forum reset init load_forum serve

load_recipes:
	python ../../manage.py project --pid test --name "Test Project" --public
	python ../../manage.py recipe --pid test --rid hello --json biostar/recipes/recipes/hello-world.hjson

load_forum:
	python manage.py populate --n_users 10 --n_posts 10 --settings ${DJANGO_SETTINGS_MODULE}

hard_reset: delete
	# Delete media and spooler.
	rm -rf ../../export/spooler/*spool*
	rm -rf ../../export/media/*

dump:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	python manage.py dumpdata --settings ${DJANGO_APP} --settings ${DJANGO_SETTINGS_MODULE} --exclude auth.permission --exclude contenttypes  > $(DUMP_FILE)
	@cp -f $(DUMP_FILE) $(BACKUP_DUMP_FILE)
	@ls -1 ../../export/database/*.json

uwsgi:
	@echo DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE}
	@echo UWSGI_INI=${UWSGI_INI}
	uwsgi --ini ${UWSGI_INI}

pg_drop:
	dropdb --if-exists ${DATABASE_NAME}
	createdb ${DATABASE_NAME}

transfer:
	python ../../manage.py migrate --settings conf.examples.pg.forum_settings
	python ../../manage.py transfer -n 300 --settings biostar.transfer.settings


config:
	(cd conf/ansible && ansible-playbook -i supportupgrade.bioconductor.org config.yml --extra-vars -v)

install:
	(cd conf/ansible && ansible-playbook -i supportupgrade.bioconductor.org install.yml --ask-become-pass --extra-vars -v)


remote_transfer:
	(cd conf/ansible && ansible-playbook -i supportupgrade.bioconductor.org transfer.yml --ask-become-pass --extra-vars -v)

deploy:
	(cd conf/ansible && ansible-playbook -i supportupgrade.bioconductor.org deploy.yml --ask-become-pass --extra-vars -v)
