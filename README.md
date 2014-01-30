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

Biostar 2.0 requires Python and a number of python packages that can be automatically installed with `pip`.
By default there are no other dependency beyond easy to install pure python packages. This means
that Biostar 2.0 will run on any platform that supports Python.

For high traffic sites third party packages are also required. These too can be automatically and easily
installed on typical Unix based systems.

Please follow the instructions in the [docs/install.md][install] file.

Deploy and customize
--------------------

Please follow the instructions in the [docs/deploy.md][deploy] file.


Migration from Biostar 1.X
--------------------------

Data migration requires a data export from the old site followed by an import into the new site.

Please follow the instructions in the [docs/migration.md][migrate] file.

[install]: docs/install.md
[migrate]: docs/migrate.md
[deploy]: docs/deploy.md
[django]: http://www.djangoproject.com/
[python]: http://www.python.org/