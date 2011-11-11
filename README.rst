BioStar Codebase
================

Introduction
-------------

BioStar codebase is a python and django based Q&A web software modeled after
the StackOverflow Q&A engine.

Our primary goal is to create a simple, flexible and extendeable 
framework.

Requirements
------------

The software requires only Python_ (2.6 or higher) to run. All other 
libraries are included in the distribution. The code will run with no
changes any operating system that supports Python.

Installation
------------

Unpack the source code archive. If you don't have django installed 
then switch to the `libs` directory and unpack the `django.zip` archive::

	$ cd libs
	$ unzip django.zip
	$ # switch back to the source directory
	$ cd ..

Quickstart
----------

From the command line execute::

    $ biostar.sh init populate run

Visit the http://localhost:8080 to view the site. Enjoy.

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

The default server will bind the all IP adapters (0.0.0.0) and port 8080. Visit http://localhost:8080 to see
interact with your version of the test server. Edit the `biostar.sh` script to override the various settings.

.. warning: The default settings will create an application with a default admin user and password!
   Modify the `main/settings.py` file to contain a different password!

Migration
---------

To migrate content from a StackExchange 1 XML datadump one needs to `import` the data. This process
may take a while. To speed up the re-import the best practice is to to `dump` the 
contents of the database it into a fixture that can be reused by 
the `populate` command. To create a new migration one would do a::

	$ biostar.sh delete init 
	$ biostar.sh import
	$ biostar.sh dump
	$ biostar.sh delete init populate

Layout
------

The Python code, templates, static content (css, images, javascript) and default 
database are found in the `main` directory. 
There is partial datadump of the existing BioStar content in the 
`import` folder. The `populate` command will load 
this data into the current database.

Other Libraries
---------------

The following software packages are being included into or being used in BioStar:

* JQuery_ for javascript and 
* `markitup`_ for rich text javascript editor. 
* `django_openid_auth`_ `python_openid`_ for openid authentication
* `pygments`_ for source code highlighting

.. _django_openid_auth: https://launchpad.net/django-openid-auth
.. _python openid: http://pypi.python.org/pypi/python-openid/
.. _pygments: http://pygments.org/

Colorscheme
-----------

  * Purple: `#8F2C47`
  * Green: `#75845C`

.. _Django: http://www.djangoproject.com/
.. _Python: http://www.python.org/
.. _JQuery: http://jquery.com/
.. _markitup: http://markitup.jaysalvat.com/home/
