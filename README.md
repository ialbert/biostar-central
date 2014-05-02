BioStar Q&A Version 2.0
=======================

BioStar is a [Python][python] and [Django][django] based Q&A software.
Our goal is to create a simple, generic, flexible and extensible Q&A
framework.

The site has been developed by scientists and for scientists. It aims to
address specific needs that scientific communities have.

This software runs several science oriented Q&A sites:

 * Biostars Bioinformatics Q&A at: https://www.biostars.org
 * Galaxy User support site: https://biostar.usegalaxy.org
 * Metabolomics Q&A: http://www.metastars.org

Features
--------

 * Post, user moderation, voting, badges, threaded discussions
 * Full email integration: import previous posts from mailing lists, support responding to posts via email

The developers of the software may be available to provide commercial level support
when deploying sites for entire organizations. Contact: admin@biostars.org

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
