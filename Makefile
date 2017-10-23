USER=www
SERVER=metabarcode.com

FASTQ_JSON = biostar/tools/fastqc/fastqc.hjson
FASTQ_TEMPL = biostar/tools/fastqc/fastqc.sh
QC_JSON =  biostar/tools/qc/qc.hjson
QC_TEMPL =  biostar/tools/qc/qc.sh
CLASSIFY_JSON = biostar/tools/classify/classify.hjson
CLASSIFY_TEMPL = biostar/tools/classify/classify.sh

EMAIL_JSON = biostar/tools/admin/send_email.hjson
EMAIL_TEMPL = biostar/tools/admin/send_email.sh


serve: init
	python manage.py runserver

uwsgi: init
	uwsgi  --ini conf/devel/devel_uwsgi.ini


delete:
	# Ensure the files that could hold secrets exist.
	touch conf/main/main_secrets.py
	touch conf/test/test_secrets.py
	# Remove the database and old media.
	rm -f export/engine.db
	rm -rf export/media/*

reset:delete init jobs adminjobs

next:
	python manage.py job --next

hello:
	python manage.py analysis --add --json biostar/tools/hello/hello4.hjson  --template biostar/tools/hello/hello4.sh --create_job
	python manage.py analysis --add --json biostar/tools/hello/hello3.hjson  --template biostar/tools/hello/hello3.sh --create_job
	python manage.py analysis --add --json biostar/tools/hello/hello2.hjson  --template biostar/tools/hello/hello2.sh --create_job
	python manage.py analysis --add --json biostar/tools/hello/hello1.hjson  --template biostar/tools/hello/hello1.sh --create_job


adminjobs:
	python manage.py analysis --add --json $(EMAIL_JSON) --template $(EMAIL_TEMPL) --analysis_usage admin  --create_job


jobs:
	python manage.py analysis --add --json $(FASTQ_JSON) --template $(FASTQ_TEMPL) --create_job
	python manage.py analysis --add --json $(QC_JSON)  --template $(QC_TEMPL) --create_job
	#python manage.py analysis --add --json $(CLASSIFY_JSON)  --template $(CLASSIFY_TEMPL) --create_job
	#python manage.py analysis --add --json $(FASTQ_JSON) --template $(FASTQ_TEMPL) --create_job
	#python manage.py analysis --add --json $(QC_JSON)  --template $(QC_TEMPL) --create_job
	#python manage.py analysis --add --json $(CLASSIFY_JSON)  --template $(CLASSIFY_TEMPL) --create_job

init:
	python manage.py collectstatic --noinput -v 0
	python manage.py migrate


test:
	python manage.py test

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

restart_nginx:
	fab -f conf/fabfile.py -H ${USER}@${SERVER} restart_nginx

restart_django:
	fab -f conf/fabfile.py -H ${USER}@${SERVER} restart_uwsgi


deploy_all: deploy_test deploy_main restart_django
