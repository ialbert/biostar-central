## Deploy

### Getting started

Thanks to the modular structure of its design Biostar is able to integrate with
a wide variety of backends and provides a number of configuration scripts and helper
methods to different deployment options.

The choices made when deploying Biostar depend on the expected levels
of traffic and number of posts that the site needs to manage. The examples that
we provide are the two extremes, some deployments may use a combination of settings from both.

Example files can be found in the `live` folder named `deploy.env` and `deploy.py`.

The basic rule is to create a settings file based on the default settings. This means that
the customized settings file will start with::

    from biostar.settings.base import *

Then subsequently override the various settings for the current deployment. For example::

    from biostar.settings.base import *
    SITE_DOMAIN = "mysite.com"
    SERVER_EMAIL = "myemail@mysite.com"

etc.

Technically a django deployment needs only a settings file, but in practice we use an environment
file to populate a shell environment and a settings file that pulls some of these variables out of
the environment.

We recommend that you start with the files in `live/deploy*` and copy them another
name. The `deploy.env` and `deploy.py` files show the minimally necessary variables
that need to be set.

    source live/deploy.env
    ./biostar.sh test

The `deploy.env` must specify the correct django settings module in this case `live.deploy` that will
load the `live/deploy.py` python module.

To run periodic scripts make sure that they load up the enviroment variables before executing the
script.

### Low traffic deployment

Suited to websites that distribute information to smaller organizations. It can be achieved
with just python based solutions. Install the dependencies with::

    pip install -r conf/requirements/deploy.txt

Copy the `live/deploy.env` and `live/deploy.py` files to a different
name/location.  For example `simple.env` and `simple.py`.
Customize these as needed. To run the site invoke the waitress server that
was installed above::

    source live/simple.env
    waitress-serve --port 8080 live.deploy.simple_wsgi:application

Create a crontab entry that updates the index every 30 minutes::

    source live/simple.env
    biostar.sh update_index

You are done.

### High traffic deployment

While not required to be turned on the site supports compressing and precompiling the site assets.
To make use of this functionality you will need to have `lessc` to be installed and you will
need to set the `USE_COMPRESSOR=True` in your settings file.

To deploy the site with `postgresql` and `elasticsearch` install the requirements::

    pip install --upgrade -r conf/requirements/deploy.txt

Start with the `conf/defaults.env` and files unde `conf/deploy/*` and customize them.
We typically copy these into the `live` folder. Rember to add an `__init__.py` file in
this folder if you want to import your settings from it.

For high performance installation we recommend deploying the production servers with
the following stack:

* Front end webserver with `nginx`
* Biostar WSGI running via `gunicorn`
* `Postgresql` as the database
* `Redis` as the job queue
* `Celery` for running the asynchronous jobs
* `Supervisord` keeping everything running
* `Elasticsearch` as the search engine

The `conf/server` folder has configuration files for `nginx`, `gunicorn` and `supervisord`.
The `conf/fabs` folder has Fabric files to automate a large number of site deployment operations.


