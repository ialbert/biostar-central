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

Migration
---------

To load content from a StackExchange 1 XML datadump one needs to *migrate* the data 
into the new schema. This is accomplished via the `main/migrate.py` script. 
This script will make use of django settings module as well. The biostar run manager's migrate 
command invokes this script. Run this script with the -h flag to see the flags it can take::

    $ python main/migrate.py -h

The result of a data migration is a json data fixture file that can be used via the *import* 
command::

    $ ./biostar.sh init import

There is an automatic account migration based on the email provided by the
OpenID provider. Only the information from a subset of well known OpenID
providers are trusted enough to allow automatic account merging. Accepted
providers are: Google, Yahoo, Myopenid, LiveJournal, Blogspot, AOL, and
Wordpress. For other users manual migration of accounts will be required.

Users listed in the Django *ADMINS* settings will have full administration privileges.

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

In Biostar there are four types of users: anonymous users, registered users, moderators and administrators.

Anonymous Users
^^^^^^^^^^^^^^^

May browse all content of a site.

Registered Users
^^^^^^^^^^^^^^^^

In addition to the privileges that anymous users have registered users  may post questions if their reputation exceeeds 
a limit (the default is zero), may post answers and comments. 

Moderator Role
^^^^^^^^^^^^^^

In addition to the privileges that registered users have moderators may edit, close and delete posts, edit user information (other than email) 
and may also suspend and reinstate users. All the actions of the moderators 
may be followed via the Moderator Log page (see About BioStar page for a link)

Administrator Role
^^^^^^^^^^^^^^^^^^

In addition to the privileges that moderators have administrators 
may promote/demote users from having moderator roles. Administrators also have 
access to the django admin interface where they may perform more database actions
than those offered via the BioStar interface..

Content Persistence
^^^^^^^^^^^^^^^^^^^

Content may be deleted (marked invisible to users) or destroyed (removed from the database).
A post submitted for deletion will be destroyed only if the author requests the deletion of
the post and the post has not collected any answers or comments. In all other cases
the post will be marked invisible to regular users.

Code Layout
-----------

The Python code, templates, static content (css, images, javascript) and default 
database are found in the *main* directory. There is partial datadump of the existing BioStar content in the 
*import* folder. The *populate* command will load 
this data into the current database.

Other Libraries
---------------

Biostar is built with open source libraries. The following software packages are used and if necessary
included with BioStar:

* Bootstrap_ as a CSS framework
* JQuery_ for javascript programming
* Less_ used for syntactically awesome css
* Coffescript_ for making javascript fun again
* markitup_ as rich text javascript editor. 
* markdown_ python library to generate HTML
* django_openid_auth_ and python_openid_ for openid authentication
* whoosh_ provides fast full text searching
* coverage_ is used to measure code coverage during testing
* prettify_ is used for syntax highlighting

.. _django_openid_auth: https://launchpad.net/django-openid-auth
.. _python_openid: http://pypi.python.org/pypi/python-openid/
.. _whoosh: https://bitbucket.org/mchaput/whoosh/wiki/Home
.. _markdown: http://www.freewisdom.org/projects/python-markdown/
.. `Python`_: http://python.org/
.. _Django: http://www.djangoproject.com/
.. _Python: http://www.python.org/
.. _JQuery: http://jquery.com/
.. _markitup: http://markitup.jaysalvat.com/home/
.. _Less: http://lesscss.org/
.. _prettify: http://code.google.com/p/google-code-prettify/
.. _Bootstrap: http://twitter.github.com/bootstrap/
