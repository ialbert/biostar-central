BioStar Q&A Version 2.0
=======================

**NOTE: This branch is under heavy development, it should not
be used by third parties until this notice disappears (posted on Jan 13, 2013)**

The version of Biostar that is currently deployed
can be found in branch `biostar1` at

https://github.com/ialbert/biostar-central/tree/biostar1

BioStar is a [Python][python] and [Django][django] based Q&A software.

Our goal is to create a simple, generic, flexible and extensible Q&A
framework.

Requirements: `Python 2.7`

Documentation
-------------

The documentation is maintained at:

http://docs.biostars.org/

The source for the documentation can be found in  the `docs` folder.

Quick Start
------------

From the biostar source directory:

    # Install the requirements.
    pip install --upgrade -r conf/requirements/base.txt

    # Load the environment variables.
    source conf/defaults.env

    # Initialize database, import test data, index for searching and run the server.
    ./biostar.sh init import index run

Visit `http://locahost:8080` to see the site loaded with default settings.

Enjoy.


[django]: http://www.djangoproject.com/
[python]: http://www.python.org/
