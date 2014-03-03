Deploy
======

Getting started
---------------

Thanks to the modular structure of its design Biostar is able to integrate with
a wide variety of backends and provides a number of configuration scripts and helper
methods to different deployment options.

The choices made when deploying Biostar depend on the expected levels
of traffic and number of posts that the site needs to manage.

The examples that we provide are the two extremes, some deployments may
use a combination of settings from both.


Low traffic deployment
-----------------------

Low concurrency: dozens of hits per minute, few emails sent.
Suited to websites that distribute information to smaller organizations.

TODO... waitress and django_static

High traffic deployment
-----------------------

TODO...

A typical deployment requires ``lessc`` to be installed and a number of other python libraries::

    pip install --upgrade -r conf/requirements/all.txt

Start with the ``conf/defaults.env`` and ``biostar/settings/deploy.py`` files and customize them.

Investigate the ``server_config`` task in ``conf/fabs/fabfile.py`` to see how we automatized the process.

There are different deployment strategies that one might follow. The sites is quite performant
and underlow concurrency the site can operate well even with the default settings of an
sqlite database running via python based webserver.

For optimal results we recommend deploying the production servers with the following stack:

* Front end webserver with ``nginx``
* Biostar WSGI running via ``gunicorn``
* ``Postgresql`` as the database
* ``Redis`` as the job queue
* ``Celery`` for running the asynchronous jobs
* ``Supervisord`` keeping everything running

The ``conf/server`` folder has configuration files for ``nginx``, ``gunicorn`` and ``supervisord``.
The ``conf/fabs`` folder has Fabric files to automate a large number of site deployment operations.


