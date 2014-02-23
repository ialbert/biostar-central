Installation
------------

Get the source and switch to the source directory

The recommended installation is via `virtualenv` and `pip`:

    pip install -r requirements/base.txt

The site manager is `biostar.sh`. This script can take one or more commands like so:

    ./biostar.sh delete init import run

Visit `http://locahost:8080` to see the site loaded with default settings.
The default admin is `foo@bar.com` password `foobar`. The default email
handler will print to the console. You can reset the password
for any user then copy paste the password reset url into the browser.

Run the manager on its own to see all the commands at your disposal:

    ./biostar.sh





To enable searching you must the content with:

    ./biostar.sh index

Deploy
------
A typical deployment requires `lessc` to be installed and a number of other python libraries.

    pip install -r conf/requirements/all.txt

Start with the `conf/defaults.env` and `biostar/settings/deploy.py` files and customize them.

Investigate the `server_config` task in `conf/fabs/fabfile.py` to see how we automatized the process.

There are different deployment strategies that one might follow. The sites is quite performant
and underlow concurrency the site can operate well even with the default settings of an
sqlite database running via python based webserver.

For optimal results we recommend deploying the production servers with the following stack:

* Front end webserver with `nginx`
* Biostar WSGI running via `gunicorn`
* `Postgresql` as the database
* `Redis` as the job queue
* `Celery` for running the asynchronous jobs
* `Supervisord` keeping everything running

The `conf/server` folder has configuration files for `nginx`, `gunicorn` and `supervisord`.
The `conf/fabs` folder has Fabric files to automate a large number of site deployment operations.

Customize
---------

See the `conf/defaults.env` for all the parameters that need to be customized.

The `SITE_STYLE_CSS` and `SITE_LOGO` settings permit loading up custom sytles. See the `/static/themes` folder
for examples.


Social Authentication
---------------------

The social logins settings will need to be initialized with the proper authentication parameters. Typically
this involves creating an application at the provider and obtaining the credentials.

See the `conf/defaults.env` for the proper variable naming.

Adding Facebook authentication:

* [Create Authentication App](http://developers.facebook.com/setup/)
* More information [Facebook Developer Resources](http://developers.facebook.com/docs/authentication/)

Adding Google authentication:

* [Google Developer Console](https://cloud.google.com/console/project)
* Create new project and copy data from credentials
* Callback must be `http://domain/accounts/google/login/callback/`

Twitter:

* Add your application at [Twitter Apps Interface](http://twitter.com/apps/)