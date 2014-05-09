.. Biostar Central documentation master file, created by
   sphinx-quickstart on Sun Mar  2 08:46:45 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Biostar: Software for building Scientific Communities
=====================================================

BioStar is a `Python <http://www.python.org/>`_ and
`Django <http://www.djangoproject.com/>`_ based Q&A software licensed under the *MIT Open Source License*.
It is a simple, generic, flexible and extensible Q&A framework.

The site has been developed by **scientists and for scientists**.
It aims to address requirements and needs that scientific communities have.

Biostar is the software that runs several science oriented Q&A sites:

 * Biostars Bioinformatics Q&A at: https://www.biostars.org
 * Galaxy User support site: https://biostar.usegalaxy.org
 * Metabolomics Q&A: http://www.metastars.org
 * Neurostars: http://www.neurostars.org


Features
--------

The site has been developed by scientists for scientists. It aims to
address specific needs that scientific communities have.

 * Standard Q&A: post questions, answers, comments, user moderation, voting, badges, threaded discussions
 * Email integration: import previous posts from mailing lists, reply to posts via email
 * RSS Planet: feed aggregation from different sources
 * External authentication: authenticate users with a different web service

Support
-------

The software is open source and free to use under the most permissible license.

The developers of the software are also available to provide commercial level support
for deploying Biostar sites for entire organizations. Contact: admin@biostars.org

Requirements: *Python 2.7*

Official documentation is located at http://docs.biostars.org

Quick Start
------------

From the biostar source directory::

    # Install the requirements.
    pip install --upgrade -r conf/requirements/base.txt

    # Load the environment variables.
    source conf/defaults.env

    # Initialize database, import test data, and run the site.
    ./biostar.sh init import index run

Visit **http://localhost:8080** to see the site loaded with default settings. Enjoy.

For more information see the documentation below:

.. toctree::
   :maxdepth: 2

   install
   manage
   deploy
   customize
   about

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

