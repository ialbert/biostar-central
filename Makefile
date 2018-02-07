USER=www
SERVER=bioinformatics.recipes

DATA_FILE=recipes-initial-data.tar.gz
DATA_DIR=/export/sites/main_data/initial
DATA_HOST=data.bioinformatics.recipes

serve: init
	python manage.py runserver


init:
	@python manage.py collectstatic --noinput -v 0
	@python manage.py migrate -v 0

uwsgi: init
	uwsgi  --ini conf/devel/devel_uwsgi.ini


install:
	pip install -r conf/python_requirements.txt
	conda config --add channels r
	conda config --add channels conda-forge
	conda config --add channels bioconda
	conda install --file conf/conda_requirements.txt -y


projects:
	python manage.py project --json initial/tutorial/tutorial-project.hjson --privacy public --jobs
	python manage.py project --root ../biostar-recipes --json projects/cookbook/cookbook-project.hjson --privacy public --jobs

	@# Create initial users
	python manage.py add_user initial/initial-users.csv
	python manage.py add_access initial/initial-access.csv

verbose:
	export DJANGO_LOG_LEVEL=DEBUG

delete:
	# Ensure the files that could hold secrets exist.
	# Remove the database and old media.
	rm -rf export/logs/*.log
	rm -rf export/spooler/*spool*
	rm -f export/database/engine.db
	rm -rf export/static/CACHE
	rm -rf export/media/*
	rm -rf *.egg
	rm -rf *.egg-info


reset: delete init projects

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

data:
	(cd export && curl http://data.bioinformatics.recipes/initial/${DATA_FILE} > ${DATA_FILE} )
	(cd export && tar zxvf ${DATA_FILE})

ftp:
	python manage.py ftp

push:
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push

