BioStar Q&A Engine
==================

**NOTE: This branch is under heavy development, it should not be used by third parties
until this notice disappears (Jan 13, 2013)**

BioStar is a [Python][python] and [Django][django] based Q&A web software.

Our primary goal is to create a simple, generic, flexible and extensible Q&A
framework.

Requirements: `Python 2.7`

Quick Start
------------

From the biostar source directory:

    pip install -r requirements/base.txt
    ./biostar.sh init import run

Visit `http://locahost:8080` to see the site loaded with default settings.

The default admin login is `foo@bar.com` with the password `foobar`.

The default email handler will print to the console.

More information on how to install, deploy and customize in [docs/install.md][install] file.

Other documentation
-------------------

* How to provide authentication from a different website. See [External authentication](docs/external.md)
* How to [migrate data from a Biostar 1.5 site](docs/migrate.md)

[django]: http://www.djangoproject.com/
[python]: http://www.python.org/