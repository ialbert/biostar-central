USER=www
SERVER=metabarcode.com

serve: init
	python manage.py runserver

install:
	pip install -r conf/python_requirements.txt
	python setup.py develop


conda:
	conda install --file conf/conda_requirements.txt -y


uwsgi:init
	uwsgi  --ini conf/devel/devel_uwsgi.ini

spool: delete init develop uwsgi

clean:
	(cd export/local && make clean)

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


reset: delete init

next:
	python manage.py job --next

hello:
	python manage.py analysis --add --json biostar/tools/hello/hello4.hjson  --template biostar/tools/hello/hello4.sh --create_job
	python manage.py analysis --add --json biostar/tools/hello/hello3.hjson  --template biostar/tools/hello/hello3.sh --create_job
	python manage.py analysis --add --json biostar/tools/hello/hello2.hjson  --template biostar/tools/hello/hello2.sh --create_job
	python manage.py analysis --add --json biostar/tools/hello/hello1.hjson  --template biostar/tools/hello/hello1.sh --create_job

jobs:
	@python manage.py analysis --add --json biostar/tools/fastqc/fastqc.hjson  --template biostar/tools/fastqc/fastqc.sh --create_job
	@python manage.py analysis --add --json biostar/tools/qc/qc.hjson  --template biostar/tools/qc/qc.sh --create_job
	@python manage.py analysis --add --json biostar/tools/classify/classify.hjson  --template biostar/tools/classify/classify.sh --create_job
	@python manage.py analysis --add --json biostar/tools/align/align.hjson  --template biostar/tools/align/align.sh --create_job
	@python manage.py analysis --add --json biostar/tools/lamar_align/lamar_align.hjson  --template biostar/tools/lamar_align/lamar_align.sh --create_job
	@python manage.py analysis --add --json biostar/tools/igv/igv.hjson  --template biostar/tools/igv/igv.sh --create_job

init:
	@python manage.py collectstatic --noinput -v 0
	@python manage.py migrate
	@python manage.py project --json initial/tutorial-project.hjson --privacy public --sticky
	@python manage.py project --json initial/cookbook-project.hjson --privacy public --sticky
	@python manage.py project --json initial/giraffe-project.hjson --sticky
	@python manage.py project --json initial/fish-project.hjson

test:
	python manage.py collectstatic --noinput -v 0
	python manage.py test -v 2 --failfast

push:
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push


test_push:test
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push

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
