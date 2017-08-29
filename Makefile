all:
	python manage.py runserver

init:
	python manage.py migrate

push:
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push
