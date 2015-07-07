[![Build Status](https://travis-ci.org/ialbert/biostar-central.svg?branch=3.0)](https://travis-ci.org/ialbert/biostar-central)

Biostar: Software for Building Scientific Communities
=====================================================

Branch 3.0 rewrite. Under heavy development. See: https://github.com/ialbert/biostar-central/issues/291

Features may not always work. Docs may be out of sync. We'll clean it up by the first release.

Quick Start
-----------

Install requirements. The site works with versions of Python 2.7 as 3.0 and above.

	pip install -r init/pip/base.txt

Run the site with default data over an `sqlite` database:

	source run/sqlite.env
	./biostar.sh delete migrate index run

To add the static pages run:

	python manage.py add_pages --dir themes/default/pages/

Visit `http://www.lvh.me:8080/` to interact with your site. It is important that even
for your localhost you use a fully qualified domain name as the group 
feature will identify groups based on the subdomain of the site.

To log in as any user note their user id.
For each user you may log in with their user id as  `email=1@foo.bar` and `password=1@foo.bar`

To run a postgresql based site use:

	source run/postgres.env
	./biostar.sh pg_create import init run

The best practice is to copy the settings into a new file for example
`mysettings.env` that are not in the repository and make use of that.
Same with the Django settings.

Documentation
-------------

* [Introduction](docs/index.md)
* [Site setup and customizations](docs/setup.md)


[django]: http://www.djangoproject.com/
[python]: http://www.python.org/
