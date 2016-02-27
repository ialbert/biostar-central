Biostar: Software for building Scientific Communities
=====================================================

[![Build Status][build-image]][build-url] 
[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)

[build-image]: https://travis-ci.org/ialbert/biostar-central.svg?branch=4.0
[build-url]: https://travis-ci.org/ialbert/biostar-central/builds

Detailed documentation: http://docs.biostars.org/

BioStar is a [Python][python] and [Django][django] based Q&A software.
It is a simple, generic, flexible and extensible Q&A framework.

The site has been developed by **scientists and for scientists**. It aims
to address the requirements and needs that scientific communities have.

Biostar is the software that runs several science oriented Q&A sites:

 * Biostars Bioinformatics Q&A at: https://www.biostars.org
 * Galaxy User support site: https://biostar.usegalaxy.org
 * Bioconductor User support site: https://support.bioconductor.org/
 * Metabolomics Q&A: http://www.metastars.org
 * Neurostars: http://www.neurostars.org


Features
--------

 * Standard Q&A: post questions, answers, comments, user moderation, voting, badges, threaded discussions
 * Email integration: import previous posts from mailing lists, reply to posts via email
 * RSS Planet: feed aggregation from different sources
 * External authentication: authenticate users with a different web service
 * Low resource utilization and easy deployment. Host an entire domain, the database, and
   all related background jobs and serve millions of page views per month on
   a computer with just 4GB of RAM.


Support
-------

The software is open source and free to use under the most permissible license.

The developers of the software are also available to provide commercial level support
for deploying Biostar sites for entire organizations. Contact: admin@biostars.org

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

Visit `http://www.lvh.me:8080` to see the site loaded with default data.
The `www.lvh.me` domain resolves to `127.0.0.0` and is your local host
with a nicer domain name. You may just as well 
use `http://localhost:8080` or `http://127.0.0.0`.

In the default site user emails are built from database ids like so :
`1@lvh.me`, `2@lvh.me`. User passwords are identical to the email. 
You may then use these to log into your test site as any of the users. 

Note that the first user always has staff level permissions and can 
access the admin interface at `http://www.lvh.me:8080/admin/`

Enjoy.


[django]: http://www.djangoproject.com/
[python]: http://www.python.org/
