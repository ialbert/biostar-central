BioStar Codebase
================

Introduction
-------------

BioStar codebase is a Python_ and Django_ based Q&A web software modeled after
the StackOverflow Q&A engine.

Our primary goal is to create a simple, generic, flexible and extendeable 
Q&A framework.

Requirements
------------

The software requires only Python_ (2.6 or higher) to run. All other 
libraries are included in the distribution. The code will run with 
no changes on any operating system that supports Python. 

The test version of this site can be seen at: http://test.biostars.org

Installation
------------

Unpack the source code archive. There are a few dependencies that
are also included with Biostar. These only need to be installed
if you don't already have them on your system.
Switch to the *libs* directory and unpack the *docutils.zip* and the *django.zip* archives::

    $ cd libs
    $ unzip docutils.zip
    $ unzip django.zip
    $ # switch back to the source directory
    $ cd ..

For faster loading performance may also want to unzip the entire `libraries.zip`
file located in the libs folder. 

Quickstart
----------

From the command line execute::

    $ ./biostar.sh init import run

Visit the http://localhost:8080 to view your site. Enjoy!

**Note** The Windows version of the biostar.sh manager has not yet
been written. The site will work just fine on Windows
but for now users will need to manually invoke the commands
present in the *biostar.sh* run manager.

Detailed Usage
--------------

There is a main run manager in the root directory::

    $ ./biostar.sh 

Execute it with no parameters for information on usage. This run manager 
can take one or more commands. For example to initialize the database then populate it with
the test data and to run the server one would invoke it in the following way::

    $ ./biostar.sh init 
    $ ./biostar.sh import
    $ ./biostar.sh run

Alternatively one may run all these commands all at once::

    $ ./biostar.sh init import run

**Note**: If database models change you must reset and reinitialize the database,
note that this will remove all existing content! The database re-initialization is
database specific, for the default sqlite deployment you can use::

    $ ./biostar.sh delete init import

The *biostar.sh* run manager to pulls in environment variables to allow you to 
customize locations/test fixtures, etc. Edit the *biostar.sh* script 
to override the various settings.

The default server will bind the all IP adapters (0.0.0.0) and port 8080. Visit http://localhost:8080 to see
interact with your version of the test server. 

Data Migration
---------------

To load content from a StackExchange 1 XML datadump one needs to *migrate* the data 
into the new schema. This is accomplished via the `migrate` command::

	$ ./biostar.sh migrate

This command in turn invokes the `main/migrate.py` script. Run this script 
(note that the Django settings need to be properly set beforehand) 
with the -h flag to see the flags it can take.::

    $ python -m main.migrate.py -h

.. note:: The `migrate` command used via the `biostar.sh` run manager makes use 
   of an in memory database as specified in the `conf/memory.env` and `conf/memory.py` files.

The result of a data migration is a compressed json data fixture file that, in turn, 
may be used via the *import* command::

    $ ./biostar.sh init import

Account migration
-----------------

There is an automatic account migration based on the email provided by the
OpenID provider. Only the information from a subset of well known OpenID
providers are trusted enough to allow automatic account merging. Accepted
providers are: Google, Yahoo, Myopenid, LiveJournal, Blogspot, AOL, and
Wordpress. For other users manual migration of accounts will be required.
Users listed in the Django *ADMINS* settings will have full administration privileges.

There is a postgresql database management script in `conf/pg-manager.sh` that is
used to facilitate data dumps and restoration.

Environment variables may be used to customize the behavior:

- `DJANGO_SETTINGS_MODULE`: the configuration module for Django
- `PYTHON`: the python executable that is to be invoked
- `FIXTURE`: output path to the (gzipped) file that will contain the data fixture
- `MIGRATE_PATH`: path to the directory that stores the StackExchange XML dump
- `MIGRATE_LIMIT`: the number of records to load from the XML dump

For a current Biostar run with about 4K users, 30K posts, 40K edits, 60K votes
generates about 300K database entries of various kinds. Data migration into a fixture
takes about 1 hour and 10Gb of RAM. This is an area that we
could do a lot better job (possibly orders of magnitude better).

The resulting data fixture is database independent and can now be loaded
into type database: sqlite, mysql, postgresql supported by Djano. For example
when loading into postgresql it takes about 2 hours and 2Gb of RAM.

Note that the databases can be dumped and restored with far fewer resources.
Exporting directly into/from postgresql for example takes less than a few
minutes.

Testing
-------

Testing also measures code coverage and therefore 
requires the coverage_ module. For your convenience this module
is included in the `libs/libraries.zip` archive. 
Install coverage_ or unzip the archive.

Testing can be initiated via the `biostar.sh` run manager::

    ./biostar.sh test

A `reports` directory will be created in the root directory
that contains html reports on the code coverage by the tests. View the `report/index.html` file.

.. _coverage: http://pypi.python.org/pypi/coverage

How the site works
-------------------

Posts may be formatted in Markdown_ (default) or ReST_ markup standards. The second format, ReST_, will be 
triggered by starting the post with the `.. rest::` directive.

User reputation is a sum of all upvotes and accepted answers that a user accumulates. Note that multiple answers
may be accepted on a question, in effect this provides the author of a question to reward twice the 
excellent answers.

In Biostar there are four types of users: anonymous users, registered users, moderators and administrators.

anonymous users
	May browse all content of a site.

registered users
	In addition to the privileges that anymous users have registered users may create new posts if their reputation 
	exceeeds a limit (the default is zero), may vote and post answers and comments. 

moderators
	In addition to the privileges that registered users have moderators may edit, close and delete posts, edit user information (other than email) 
	and may also suspend and reinstate users. All the actions of the moderators 
	may be followed via the Moderator Log page (see About BioStar page for a link)

administrators
	In addition to the privileges that moderators have administrators 
	may promote/demote users from having moderator roles. Administrators also have 
	access to the django admin interface where they may perform more database actions
	than those offered via the BioStar interface..

Content Persistence
^^^^^^^^^^^^^^^^^^^

Content may be deleted (marked invisible to users) or destroyed (removed from the database).

A post submitted for deletion will be destroyed only if the author requests the deletion 
and the post does not have any followups (answers/comments) associated with it. Deleted top level posts 
are marked invisible to regular users.

Code Layout
-----------

The Python code, templates, static content (css, images, javascript) and default 
database are found in the *main* directory. There is partial datadump of the existing BioStar content in the 
*import* folder. The *import* command will load this data into the current database.

Other Libraries
---------------

Biostar is built with open source libraries. The following software packages are used and 
if necessary included and distributed with BioStar:

* Bootstrap_ as a CSS framework
* JQuery_ for javascript programming
* Less_ used for syntactically awesome css
* markitup_ as rich text javascript editor. 
* python-markdown_ python library to convert Markdown_ to  HTML
* docutils_ is used to convert ReST_ to HTML
* django_openid_auth_ and python_openid_ for openid authentication
* whoosh_ provides fast full text searching
* coverage_ is used to measure code coverage during testing
* prettify_ is used for syntax highlighting


.. _django_openid_auth: https://launchpad.net/django-openid-auth
.. _python_openid: http://pypi.python.org/pypi/python-openid/
.. _whoosh: https://bitbucket.org/mchaput/whoosh/wiki/Home
.. _python-markdown: http://www.freewisdom.org/projects/python-markdown/
.. `Python`_: http://python.org/
.. _Django: http://www.djangoproject.com/
.. _Python: http://www.python.org/
.. _JQuery: http://jquery.com/
.. _markitup: http://markitup.jaysalvat.com/home/
.. _Less: http://lesscss.org/
.. _prettify: http://code.google.com/p/google-code-prettify/
.. _Bootstrap: http://twitter.github.com/bootstrap/
.. _docutils: http://docutils.sourceforge.net/docs/user/rst/quickstart.html
.. _ReST: http://docutils.sourceforge.net/docs/user/rst/quickstart.html
.. _Markdown: http://en.wikipedia.org/wiki/Markdown