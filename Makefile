USER=www
SERVER=metabarcode.com

serve: init
	python manage.py runserver

uwsgi: init
	uwsgi  --ini live/devel_uwsgi.ini


delete:
	rm -f export/engine.db
	rm -rf media/*

reset:delete init


init:
	python manage.py collectstatic --noinput -v 0
	python manage.py migrate

test:
	python manage.py test

push:
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push

testdata:
	mkdir -p tmp/refs/
	curl  http://iris.bx.psu.edu/projects/metabarcode-data/data.tar.gz > tmp/data.tar.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/sampleinfo.txt > tmp/sampleinfo.txt
	curl http://iris.bx.psu.edu/projects/metabarcode-data/data/1-SarriPal_S9_L001_R1_001.fastq.gz > tmp/1-SarriPal_S9_L001_R1_001.fastq.gz
	curl http://iris.bx.psu.edu/projects/metabarcode-data/centrifuge_fish_idx.tar.gz >tmp/refs/centrifuge_fish_idx.tar.gz
	(cd tmp/refs  && tar -xzvf centrifuge_fish_idx.tar.gz) 

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
