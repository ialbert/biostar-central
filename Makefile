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

serve_pg:
	python manage.py runserver --settings conf.postgres.postgres_settings

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

tutorial:
	# Perform a full delete of all content.
	python manage.py data --pid test --path /export/data/usfish/runs/test/ --name "Test Sequences" --type FASTQ
	python manage.py data --pid test --path /export/data/usfish/reference/fish-accession.fa --type FASTA


full_delete: delete
	# Perform a full delete of all content.
	rm -rf export/spooler/*spool*
	rm -rf export/media/*

# Resets the site without removing jobs.
reset: delete init
    # Initializes the test project.
	python manage.py project --pid test --name "Test Project" --public
	python manage.py recipe --pid test --rid hello --json biostar/engine/recipes/hello-world.hjson
	python manage.py recipe --pid test --rid kraken2 --json ~/book/biostar-handbook-2/recipes/work/classify/kraken2.hjson

# Resets site and loads existing fixture
reset_load: delete init loaddata

# Removes all content from the site.
hard_reset: full_delete init

loaddata:
	python manage.py loaddata --ignorenonexistent $(DUMP_FILE)

dumpdata:
	python manage.py dumpdata --exclude auth.permission --exclude contenttypes --exclude transfer > $(DUMP_FILE)
	cp -f $(DUMP_FILE) $(BACKUP_DUMP_FILE)
	# Produce a datadump count as a reminder.
	@ls -1 export/database/*.json | wc -l

uwsgi:
	uwsgi  --ini conf/devel/devel_uwsgi.ini

serve_forum:
	python manage.py runserver --settings=biostar.forum.settings

pgserve_forum:
	python manage.py runserver --settings=biostar.forum.postgres_settings


projects:
	# Load projects from json files in --dir
	python manage.py api create --dir ../biostar-recipes/projects

	# Load data associated with projects
	python manage.py api create --json ../biostar-recipes/data.hjson --is_data


recipes2: projects

	# Load recipes from json files in --dir, start jobs on some
	python manage.py api create --dir ../biostar-recipes/recipes/cookbook/ --create_jobs
	python manage.py api create --dir ../biostar-recipes/recipes/tutorial --create_jobs
	python manage.py api create --dir ../biostar-recipes/recipes/fish/
	python manage.py api create --dir ../biostar-recipes/recipes/giraffe/
	python manage.py api create --dir ../biostar-recipes/recipes/handbook/
	python manage.py api create --dir ../biostar-recipes/recipes/mothur/
	python manage.py api create --dir ../biostar-recipes/recipes/trout/

	# Create initial users
	python manage.py add_user --fname initial/initial-users.csv
	python manage.py add_access initial/initial-access.csv

# Load all recipes
recipes: tutorial
	python manage.py project --root ../biostar-recipes --json projects/cookbook-project.hjson --privacy public --jobs
	python manage.py project --root ../biostar-recipes --json projects/mothur-project.hjson --privacy public
	python manage.py project --root ../biostar-recipes --json projects/giraffe-project.hjson --privacy public
	python manage.py project --root ../biostar-recipes --json projects/handbook-project.hjson --privacy public
	python manage.py project --root ../biostar-recipes --json projects/usfish-project.hjson --privacy public
	python manage.py project --root ../biostar-recipes --json projects/trout-project.hjson

	# Create initial users
	python manage.py add_user --fname initial/initial-users.csv
	python manage.py add_access initial/initial-access.csv

verbose:
	# Makes logging more verbose.
	export DJANGO_LOG_LEVEL=DEBUG

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
	python manage.py load --root initial/export-1000 --users users.txt --n 1000 || true

	python manage.py load --root initial/export-1000 --posts posts.txt  --n 1000 || true

	python manage.py load --root initial/export-1000 --votes votes.txt  --n 1000 || true


deploy_psu:
	(cd conf/ansible && ansible-playbook -i hosts-psu server_deploy.yml --ask-become-pass --extra-vars)

deploy_www:
	(cd conf/ansible && ansible-playbook -i hosts server_deploy.yml --ask-become-pass --extra-vars -v)

