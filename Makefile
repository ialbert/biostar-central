# Data dump files.
DUMP_FILE=export/database/db.json

# Backup file.
BACKUP_DUMP_FILE=export/database/db.backup.`date +'%Y-%m-%d-%H%M'`.json

# Default settings module.
DJANGO_SETTING_MODULE := biostar.engine.settings

# Default app.
DJANGO_APP := biostar.engine

# Database name
DATABASE_NAME := database.db

all: engine serve

accounts:
	$(eval DJANGO_SETTING_MODULE := biostar.accounts.settings)
	$(eval DJANGO_APP := biostar.accounts)

	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}

emailer:
	$(eval DJANGO_SETTING_MODULE := biostar.emailer.settings)
	$(eval DJANGO_APP := biostar.emailer)

	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}

pg:
	$(eval DJANGO_SETTING_MODULE := conf.examples.postgres.postgres_settings)
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}

message:
	$(eval DJANGO_SETTING_MODULE := biostar.message.settings)
	$(eval DJANGO_APP := biostar.message)

	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}

engine:
	$(eval DJANGO_SETTING_MODULE := biostar.engine.settings)
	$(eval DJANGO_APP := biostar.engine)
	$(eval UWSGI_INI := conf/uwsgi/engine_uwsgi.ini)

	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}

forum:
	$(eval DJANGO_SETTING_MODULE := biostar.forum.settings)
	$(eval DJANGO_APP := biostar.forum)
	$(eval UWSGI_INI := conf/uwsgi/forum_uwsgi.ini)

	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	@echo DATABASE_NAME=${DATABASE_NAME}

serve:
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	python manage.py runserver --settings ${DJANGO_SETTING_MODULE}

init:
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	python manage.py collectstatic --noinput -v 0  --settings ${DJANGO_SETTING_MODULE}
	python manage.py migrate -v 0  --settings ${DJANGO_SETTING_MODULE}

load:
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	python manage.py loaddata --ignorenonexistent --settings ${DJANGO_SETTING_MODULE} $(DUMP_FILE)

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

test:
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	python manage.py test ${DJANGO_APP} --settings ${DJANGO_SETTING_MODULE} -v 2 --failfast
	coverage run manage.py test ${DJANGO_APP} --settings ${DJANGO_SETTING_MODULE} -v 2 --failfast
	coverage html --skip-covered

test_all:
	python manage.py test --settings biostar.test.test_settings -v 2 --failfast
	coverage run manage.py test --settings biostar.test.test_settings -v 2 --failfast
	coverage html --skip-covered

projects:
	python manage.py project --pid test --name "Test Project" --public
	python manage.py recipe --pid test --rid hello --json biostar/engine/recipes/hello-world.hjson

populate:
	python manage.py populate --settings ${DJANGO_SETTING_MODULE}

hard_reset: delete
	# Delete media and spooler.
	rm -rf export/spooler/*spool*
	rm -rf export/media/*

dump:
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	python manage.py dumpdata --settings ${DJANGO_APP} --settings ${DJANGO_SETTING_MODULE} --exclude auth.permission --exclude contenttypes  > $(DUMP_FILE)
	@cp -f $(DUMP_FILE) $(BACKUP_DUMP_FILE)
	@ls -1 export/database/*.json

uwsgi:
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo UWSGI_INI=${UWSGI_INI}
	uwsgi --ini ${UWSGI_INI}

drop_create:
	dropdb --if-exists ${DATABASE_NAME}
	createdb ${DATABASE_NAME}

transfer:
	python manage.py migrate --settings conf.examples.postgres.transfer_settings
	python manage.py transfer -n 300 --settings conf.examples.postgres.transfer_settings

next:
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	python manage.py job --next --settings ${DJANGO_SETTING_MODULE}

push:
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push

deploy:
	(cd conf/ansible && ansible-playbook -i hosts server_deploy.yml --ask-become-pass --extra-vars -v)
