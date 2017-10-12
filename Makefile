USER=www
SERVER=metabarcode.com

serve: init
	python manage.py runserver

reset:
	rm -f export/engine.db
	rm -rf media/*

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
