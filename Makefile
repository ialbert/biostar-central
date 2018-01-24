USER=www
SERVER=bioinformatics.recipes

DATA_FILE=recipes-initial-data.tar.gz
DATA_DIR=/export/sites/main_data/initial
DATA_HOST=data.bioinformatics.recipes

serve: init
	python manage.py runserver

install:
	pip install -r conf/python_requirements.txt
	python setup.py develop

conda:
	conda config --add channels r
	conda config --add channels conda-forge
	conda config --add channels bioconda
	conda install --file conf/conda_requirements.txt -y


verbose:
	export DJANGO_LOG_LEVEL=DEBUG

uwsgi:init
	uwsgi  --ini conf/devel/devel_uwsgi.ini

spool: delete init develop uwsgi

clean:
	(cd export/local && make clean)

data_pack:
	(cd export/local && make  data_pack)

data_push:
	(cd export/local && make  data_push)

data_pull:
	(cd export/local && make  data_pull)

testdata:
	(cd export/local && make  all)

delete:
	# Ensure the files that could hold secrets exist.
	touch conf/main/main_secrets.py
	touch conf/test/test_secrets.py
	# Remove the database and old media.
	rm -rf export/logs/*.log
	rm -rf export/spooler/*spool*
	rm -f export/database/engine.db
	rm -rf export/static/CACHE
	rm -rf export/media/*
	rm -rf *.egg
	rm -rf *.egg-info


postgres:
	#dropdb --if-exists testbuddy_engine
	#createdb testbuddy_engine
	#python manage.py migrate
	python manage.py test --settings conf.postgres.postgres_settings --failfast


#reset: delete init tutorial cookbook fish giraffe mothur users

reset: delete init tutorial cookbook giraffe users


next:
	python manage.py job --next

users:
	@# Create the initial users.
	@python manage.py add_user initial/initial-users.csv

	@# Create initial access of users.
	@python manage.py add_access initial/initial-access.csv

init:
	@python manage.py collectstatic --noinput -v 0
	@python manage.py migrate -v 0

# Initializes the tutorial projects.
tutorial:
	python manage.py project --json initial/tutorial/tutorial-project.hjson --privacy public --jobs

cookbook:
	#python manage.py project --root ../biostar-recipes --json projects/cookbook/cookbook-project.hjson --privacy public --jobs

	#@python manage.py project --json initial/cookbook-project.hjson --privacy public --sticky
	#@python manage.py project --json initial/biostar-handbook.hjson --privacy public --sticky

fish:
	python manage.py project --json initial/fish-project.hjson

giraffe:
	python manage.py project --root ../biostar-recipes --json projects/giraffe/giraffe-project.hjson --sticky --privacy public

mothur:
	python manage.py project --root ../biostar-recipes --json projects/metagenome/mothur-project.hjson --privacy public

test:
	python manage.py collectstatic --noinput -v 0
	python manage.py test -v 2 --failfast
	coverage run manage.py test
	coverage html --skip-covered

data:
	(cd export && curl http://data.bioinformatics.recipes/initial/${DATA_FILE} > ${DATA_FILE} )
	(cd export && tar zxvf ${DATA_FILE})

ftp:
	python biostar/ftpserver/server.py

push:
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push

test_push:test push

deploy_test:
	fab -f conf/fabfile.py -H ${USER}@${SERVER} deploy_test

deploy_main:
	fab -f conf/fabfile.py -H ${USER}@${SERVER} deploy_main

deploy_psu:
	fab -f conf/fabfile.py -H www@psu.bioinformatics.recipes deploy_psu

restart_nginx:
	fab -f conf/fabfile.py -H ${USER}@${SERVER} restart_nginx

restart_django:
	fab -f conf/fabfile.py -H ${USER}@${SERVER} restart_uwsgi

deploy_all: deploy_test deploy_main restart_django
