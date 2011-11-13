BioStar Codebase
================

Introduction
-------------

BioStar codebase is a Python_ and Django_ based Q&A web software modeled after
the StackOverflow Q&A engine.

Our primary goal is to create a simple, generic, flexible and extendeable 
Q&A framework.

A demo version of this site can be seen at: http://www.biostars.org

**Note**: Please note that software is not yet robust enough for public
deployment! 

Requirements
------------

The software requires only Python_ (2.6 or higher) to run. All other 
libraries are included in the distribution. The code will run with 
no changes on all operating system that support Python. 

Installation
------------

Unpack the source code archive. If you don't have django installed 
then switch to the *libs* directory and unpack the *django.zip* archive (included
for convenience)::

	$ cd libs
	$ unzip django.zip
	$ # switch back to the source directory
	$ cd ..

Quickstart
----------

From the command line execute::

    $ biostar.sh init populate run

Visit the http://localhost:8080 to view the site. Enjoy.

**Note** The Windows version of the biostar.sh manager has not yet
been written. The site will work just fine on Windows
but for now users will need to manually invoke the commands
present in the *biostar.sh* run manager (TODO).

Detailed Usage
--------------

There is a main run manager in the root directory::

    $ biostar.sh 

Execute it with no parameters for information on usage. This run manager 
can take one or more commands. For example to initialize the database then populate it with
the test data and to run the server one would invoke it in the following way::

    $ biostar.sh init 
    $ biostar.sh populate
    $ biostar.sh run

Alternatively one may run all these commands all at once::

    $ biostar.sh init populate server

**Note**: If the database models change you must reset and reinitialize your database::

    $ biostar.sh delete init populate

The *biostar.sh* run manager to pulls in environment variables to allow you to 
customize locations/test fixtures, etc. View the source for well documented examples.

The default server will bind the all IP adapters (0.0.0.0) and port 8080. Visit http://localhost:8080 to see
interact with your version of the test server. Edit the *biostar.sh* script to override the various settings.

**Warning**: The default settings will create an application with a default admin user and password!
Modify the *main/settings.py* file to contain a different password!

Migration
---------

To migrate content from a StackExchange 1 XML datadump one needs to *import* the data. This process
may take a while (TODO: speeding up import via bulk loading).
To facilitate the re-import the best practice is that after an *import* takes
place site admins should *dump* the contents of the database it into a fixture
that can be reused by the *populate* command. The data fixture loading is
substantially faster. Therefore to create a new migration one would do the
following::

	$ biostar.sh import
	$ biostar.sh delete init populate

The *ALLOW_MIGRATION* setting will enable a single automatic account migration
based on the email provided by the OpenID provider. Only the information
from a subset of well known OpenID providers are trusted enough
to automatically merge accounts. For other users manual migration of accounts
will be required.

Users in *ADMINS* settings will automatically obtain full administration privileges and
may log into the *admin* site using the *SECRET_KEY* as their password.

Layout
------

The Python code, templates, static content (css, images, javascript) and default 
database are found in the *main* directory. 
There is partial datadump of the existing BioStar content in the 
*import* folder. The *populate* command will load 
this data into the current database.

Other Libraries
---------------

The following software packages are being included into or being used in BioStar:

* JQuery_ for javascript and 
* `markitup`_ for rich text javascript editor. 
* `django_openid_auth`_ and `python_openid`_ for openid authentication
* `pygments`_ for source code highlighting

.. _django_openid_auth: https://launchpad.net/django-openid-auth
.. _python_openid: http://pypi.python.org/pypi/python-openid/
.. _pygments: http://pygments.org/

Colorscheme
-----------

  * Purple: `#8F2C47`
  * Green: `#75845C`

.. _Django: http://www.djangoproject.com/
.. _Python: http://www.python.org/
.. _JQuery: http://jquery.com/
.. _markitup: http://markitup.jaysalvat.com/home/
