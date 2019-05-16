# Data dump files.
DUMP_FILE=export/database/db.json

# Backup file.
BACKUP_DUMP_FILE=export/database/db.backup.`date +'%Y-%m-%d-%H%M'`.json

# Default settings module.
DJANGO_SETTING_MODULE := biostar.engine.settings

# Default app.
DJANGO_APP := biostar.engine

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
	$(eval DJANGO_SETTING_MODULE := conf.postgres.postgres_settings)
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}

message:
	$(eval DJANGO_SETTING_MODULE := biostar.message.settings)
	$(eval DJANGO_APP := biostar.message)

	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}

engine:
	$(eval DJANGO_SETTING_MODULE := biostar.engine.settings)
	$(eval DJANGO_APP := biostar.engine)

	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}

forum:
	$(eval DJANGO_SETTING_MODULE := biostar.forum.settings)
	$(eval DJANGO_APP := biostar.forum)

	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}

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
	rm -f export/database/engine.db
	rm -rf export/static/CACHE
	rm -rf *.egg
	rm -rf *.egg-info

# Resets the site without removing jobs.
reset: delete init
    # Initializes the test project.

test:
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	@echo DJANGO_APP=${DJANGO_APP}
	python manage.py collectstatic --noinput -v 0 --settings ${DJANGO_SETTING_MODULE}
	python manage.py test ${DJANGO_APP} --settings ${DJANGO_SETTING_MODULE} -v 2 --failfast
	coverage run manage.py test ${DJANGO_APP} --settings ${DJANGO_SETTING_MODULE} -v 2 --failfast
	coverage html --skip-covered

projects:
	python manage.py project --pid test --name "Test Project" --public
	python manage.py recipe --pid test --rid hello --json biostar/engine/recipes/hello-world.hjson

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
	uwsgi --ini conf/devel/devel_uwsgi.ini


pg_forum:
	python manage.py runserver --settings=biostar.forum.postgres_settings

pg_drop:
	dropdb --if-exists engine.db

pg_create:
	#dropdb --if-exists engine.db
	createdb engine.db
	python manage.py migrate --settings conf.postgres.postgres_settings
	python manage.py test --settings conf.postgres.postgres_settings --failfast

postgres:
	python manage.py migrate --settings conf.postgres.postgres_settings
	python manage.py test --settings conf.postgres.postgres_settings --failfast

next:
	@echo DJANGO_SETTING_MODULE=${DJANGO_SETTING_MODULE}
	python manage.py job --next --settings ${DJANGO_SETTING_MODULE}

push:
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push

deploy:
	(cd conf/ansible && ansible-playbook -i hosts server_deploy.yml --ask-become-pass --extra-vars -v)

