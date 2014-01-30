BioStar Q&A Engine
==================

**NOTE: This branch is under heavy development, it should not be used by third parties
until this notice disappears (Jan 13, 2013)**

BioStar is a [Python][python] and [Django][django] based Q&A web software.

Our primary goal is to create a simple, generic, flexible and extensible Q&A
framework.

Requirements: `Python 2.7`

Installation
------------

The recommended installation is via `virtualenv`.

	pip install -r requirements/base.txt


Migration from Biostar 1.X
--------------------------

Data migration requires a data export from the old site followed by an import into the new site.

Please follow the instructions in the [docs/migration.md][migration] file.

[migration]: tree/master/docs/migration.md
[django]: http://www.djangoproject.com/
[python]: http://www.python.org/