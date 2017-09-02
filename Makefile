USER=www
SERVER=metabarcode.com

serve:
	python manage.py runserver

reset:
	rm -f export/engine.db

init:
	python manage.py collectstatic --noinput
	python manage.py migrate

push:
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push

test_pull:
	fab -f conf/fabfile.py -H ${USER}@${SERVER} test_pull

main_pull:
	fab -f conf/fabfile.py -H ${USER}@${SERVER} main_pull

restart_nginx:
	fab -f conf/fabfile.py -H ${USER}@${SERVER} restart_nginx

restart_uwsgi:
	fab -f conf/fabfile.py -H ${USER}@${SERVER} restart_uwsgi