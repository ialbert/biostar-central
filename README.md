## Biostar: Building Scientific Communities

[![Build Status][build-image]][build-url] 
[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)

[build-image]: https://travis-ci.org/ialbert/biostar-central.svg?branch=4.0
[build-url]: https://travis-ci.org/ialbert/biostar-central/builds

BioStar is a [Python][python] and [Django][django] based Q&A software.
It is a simple, generic, flexible and extensible Q&A framework.

The site has been developed by **scientists and for scientists**. It aims
to address the requirements and needs that scientific communities have.
Biostar is used to run several science oriented Q&A sites:

 * Biostars Bioinformatics Q&A at: https://www.biostars.org
 * Galaxy User Support: https://biostar.usegalaxy.org
 * Bioconductor User Support: https://support.bioconductor.org/
 * Neuroinformatics Q&A: http://www.neurostars.org

The software is open source and free to use under the MIT License.

### Features

* Q&A: questions, answers, comments, user moderation, voting, reputation, badges, threaded discussions
* RSS Planet: feed aggregation from different sources
* External authentication: authenticate users with a different web service
* Email integration: import previous posts from mailing lists 
* Low resource utilization and easy deployment. 

Requirements: `Python 2.7`

### Documentation

The documentation:

* [Install](docs/install.md)
* [Manage](docs/manage.md)
* [Deploy](docs/deploy.md)
* [Customize](org/bioconductor/README.md)

The source for the documentation can be found in  the [docs](./docs) folder.

### Quick Start

From the biostar source directory:

    # Install the requirements.
    pip install --upgrade -r conf/requirements/base.txt

    # Initialize database, import test data, index for searching and run the server.
    ./biostar.sh init import index run

Visit `http://www.lvh.me:8080` to see the site loaded with demo data.

The `www.lvh.me` domain resolves to `127.0.0.1` your local host 
with a proper domain name. You may just as well use `http://localhost:8080` or `http://127.0.0.1`.

In the demo site the user emails are built from the database ids like so: `1@lvh.me`, `2@lvh.me`.
The demo user passwords are identical to the emails
and you may use these to log into your test site as any of the users.

The user with the email `1@lvh.me` has staff level permissions and
can also access the admin interface at `http://www.lvh.me:8080/admin/`.

Enjoy.

### Development

Biostar versions and upgrade path: https://github.com/ialbert/biostar-central/issues/400

### Support

We may be able to provide support for organizations or institutions. 
For more information contact **admin@biostars.org**

[django]: http://www.djangoproject.com/
[python]: http://www.python.org/

### Citing

* Parnell LD, Lindenbaum P, Shameer K, Dall'Olio GM, Swan DC, et al.
  [2011 BioStar: An Online Question & Answer Resource for the Bioinformatics Community.] (http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002216)
  PLoS Comput Biol 7(10): e1002216. doi:10.1371/journal.pcbi.1002216

### Contributors

List of contributors: https://github.com/ialbert/biostar-central/graphs/contributors