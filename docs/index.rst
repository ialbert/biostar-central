.. Biostar Central documentation master file, created by
   sphinx-quickstart on Sun Mar  2 08:46:45 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Biostar Central's documentation!
===========================================

BioStar is a `Python <http://www.python.org/>`_ and
`Django <http://www.djangoproject.com/>`_ based Q&A software.

Our goal is to create a simple, generic, flexible and extensible Q&A
framework.

Requirements: `Python 2.7`

Quick Start
------------

From the biostar source directory::

    # Install the requirements.
    pip install --upgrade -r conf/requirements/base.txt

    # Load the environment variables.
    source conf/defaults.env

    # Initialize database, import test data, and run the site.
    ./biostar.sh init import run

Visit ``http://locahost:8080`` to see the site loaded with default settings.

For more information see the documentation.

.. toctree::
   :maxdepth: 2

   install
   external
   migrate

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

