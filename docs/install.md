Installation
------------

Get the source and switch to the source directory

The recommended installation is via `virtualenv` and `pip`:

    pip install -r requirements/base.txt

The site manager is `biostar.sh`. This command can take one or more commands like so:

    ./biostar.sh init import run

Visit `http://locahost:8080` to see the site loaded with default settings.

To enable searching you must the content with:

    ./biostar.sh index

The next step is to [deploy Biostar][deploy].

[deploy]: docs/deploy.md
