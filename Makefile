# Test data specific variables.
DATA_FILE=recipes-initial-data.tar.gz
DATA_DIR=/export/sites/main_data/initial
DATA_HOST=data.bioinformatics.recipes

# Data dump files.
DUMP_FILE=export/database/engine.json
BACKUP_DUMP_FILE=export/database/engine_`date +'%s'`.json

all: serve

serve:
	python manage.py runserver

init:
	python manage.py collectstatic --noinput -v 0
	python manage.py migrate -v 0

delete:
	# Delete the database,logs and CACHE files
	# Keeps media and spool
	rm -rf export/logs/*.log
	rm -f export/database/engine.db
	rm -rf export/static/CACHE
	rm -rf *.egg
	rm -rf *.egg-info

full_delete: delete
	# Perform a full delete of all content.
	rm -rf export/spooler/*spool*
	rm -rf export/media/*

# Resets the site without removing jobs.
reset: delete init loaddata

# Removes all content from the site.
hard_reset: full_delete init

loaddata:
	python manage.py loaddata $(DUMP_FILE)

backup:
	python manage.py dumpdata --exclude auth.permission --exclude contenttypes > $(DUMP_FILE)
	cp -f $(DUMP_FILE) $(BACKUP_DUMP_FILE)
	# Produce a datadump count as a reminder.
	@ls -1 export/database/*.json | wc -l

dumpdata: backup

uwsgi:
	uwsgi  --ini conf/devel/devel_uwsgi.ini

serve_forum:
	python manage.py runserver --settings=biostar.forum.settings

install:
	pip install -r conf/python_requirements.txt
	conda config --add channels r
	conda config --add channels conda-forge
	conda config --add channels bioconda
	conda install --file conf/conda_requirements.txt -y

tutorial:
	python manage.py project --root ../biostar-recipes --json projects/tutorial-project.hjson --privacy public --jobs

# Load all recipes
recipes: tutorial
	python manage.py project --root ../biostar-recipes --json projects/cookbook-project.hjson --privacy public --jobs
	python manage.py project --root ../biostar-recipes --json projects/mothur-project.hjson --privacy public
	python manage.py project --root ../biostar-recipes --json projects/giraffe-project.hjson --privacy public
	python manage.py project --root ../biostar-recipes --json projects/handbook-project.hjson --privacy public
	python manage.py project --root ../biostar-recipes --json projects/usfish-project.hjson --privacy public
	python manage.py project --root ../biostar-recipes --json projects/trout-project.hjson

	# Create initial users
	python manage.py add_user initial/initial-users.csv
	python manage.py add_access initial/initial-access.csv

verbose:
	# Makes logging more verbose.
	export DJANGO_LOG_LEVEL=DEBUG

postgres:
	#dropdb --if-exists testbuddy_engine
	#createdb testbuddy_engine
	#python manage.py migrate
	python manage.py test --settings conf.postgres.postgres_settings --failfast

next:
	python manage.py job --next

test:
	python manage.py collectstatic --noinput -v 0
	python manage.py test -v 2 --failfast
	coverage run manage.py test
	coverage html --skip-covered

ftp:
	python manage.py ftp

push:
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push

# Biostar data load
biostar_load:

	@tar -xvzf initial/initial-posts.tar.gz --directory initial

	@# Load initial users first
	python manage.py load --root initial/export-100 --users users.txt --n 100 || true

	python manage.py load --root initial/export-100 --posts posts.txt  --n 100 || true

	python manage.py load --root initial/export-100 --votes votes.txt  --n 100 || true


deploy_psu:
	(cd conf/ansible && ansible-playbook -i hosts-psu server_deploy.yml --ask-become-pass --extra-vars)

deploy_www:
	(cd conf/ansible && ansible-playbook -i hosts server_deploy.yml --ask-become-pass --extra-vars -v)

