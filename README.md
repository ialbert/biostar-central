[![Build Status](https://travis-ci.org/ialbert/biostar-central.svg?branch=3.0)](https://travis-ci.org/ialbert/biostar-central)

Biostar: Software for Building Scientific Communities
=====================================================

Branch 3.0 rewrite, under heavy development. See: https://github.com/ialbert/biostar-central/issues/291

Quick Start
-----------

Install requirements. The site needs [Python 2.7][python] installed:

	pip install -r conf/pip/base.txt

Run the site with default data over an `sqlite` database:

	source run/sqlite.env
	./biostar.sh delete import init index run

Visit http://localhost:8080 to interact with your site. To log in as any user note their user id.
For each user you may log in with their user id as  `email=1@foo.bar` and `password=1@foo.bar`

To run a postgresql based site use:

	source run/postgres.env
	./biostar.sh pg_create import init run

The best practice is to copy the settings into a new file for example
`mysettings.env` that are not in the repository and make use of that.
Same with the Django settings. See more in the documentation.


[django]: http://www.djangoproject.com/
[python]: http://www.python.org/
