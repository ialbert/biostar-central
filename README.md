BioStar Codebase
================

Introduction
-------------

BioStar is a [Python][python] and [Django][django] based Q&A web software
modeled after the StackOverflow Q&A engine.

Our primary goal is to create a simple, generic, flexible and extendeable Q&A
framework. 

Requirements
------------

The software requires only [Python][python] (version 2.7) to run. Most of the
other libraries except `django-celery` are included in the distribution. This latter
needs to be present on the installed computer.

The code will run with no changes
on any operating system that supports Python. The software runs the **BioStar
Bioinformatics Q&A** site at: http://www.biostars.org

Installation
------------

Get the source from github

    $ git clone https://github.com/ialbert/biostar-central.git

Setup a virtual environment (if not installing the site system-wide):

    $ curl -s https://raw.github.com/brainsik/virtualenv-burrito/master/virtualenv-burrito.sh | $SHELL

Install dependencies:

    $ pip install -r requirements.txt

Quickstart
----------

You will need to initialize a number of environment variables then run the
server:

    $ source conf/default.env
    $ ./biostar.sh init import run

Visit the http://localhost:8080 to view your site. Enjoy!

**Note** The Windows version of the biostar.sh manager has not yet been
written. The site will work just fine on Windows but for now users will need to
manually invoke the commands present in the *biostar.sh* run manager.

Usage
-----

In general you will first need to initialize environment variables then start
the server. There are a number of example environments in the `conf` folder.
The default environment can be initialized via:

    $ source conf/default.env

There is a run manager in the root directory:

    $ ./biostar.sh 

Execute it with no parameters for information on usage. This run manager can
take one or more commands. For example to initialize the database then populate
it with the test data and to run the server one would invoke it in the
following way:

    $ ./biostar.sh init 
    $ ./biostar.sh import
    $ ./biostar.sh run

Alternatively one may run all these commands all at once:

    $ ./biostar.sh init import run

**Note**: If database models change you must reset and reinitialize the
database, note that this will remove all existing content! The database
re-initialization is database specific, for the default sqlite deployment you
can use:

    $ ./biostar.sh delete init import run

The *biostar.sh* run manager to pulls in environment variables to allow you to
customize locations/test fixtures, etc. Edit the *biostar.sh* script to
override the various settings.

Search requires indexing that is disabled during a migration or import. To
enable the search you will need to manually trigger the indexing via::

    $ ./biostar.sh index

The default server will bind the all IP adapters (0.0.0.0) and port 8080. Visit
http://localhost:8080 to see interact with your version of the test server. 

There are commands to support the Postgresql database. The commands are:

    pgdrop pgdump pgreset pgimport <your-sql-file>

Most operations are customized via environment variables. To show their current
settings use::

    $ ./biostar.sh env

Settings
--------

Custom settings should import all settings as `from main.settings import *`
then override those values that need to change. For an example see the
`conf/demo.env` and `conf/demo.py` files.  A number of helper methods can be
found in the `main/bin` directory, when the python import paths are set
properly these can be invoked as:

  - `python -m main.bin.postadd` to import posts into the main site
  - `python -m main.bin.postmod` modifies a post in the main site
  - `python -m main.bin.useradd` to add users
  - `python -m main.bin.usermod` to edit users
  - `python -m main.bin.extract` extracts the value of a settings parameter
  - `python -m main.bin.patch` is used to change various values after migration
  - `python -m main.bin.sitemap` generates a site map into the export directory

The best practice is to establish separate environments and settings files.
The `python -m main.bin.extract` command can be used to fill an environment
variable from a Django settings file. See the `conf/demo.env` file. Typically
use is to `source conf/default.env` then source a second smaller file that
overrides just a few parameters.

To add a new user and a post by this user you can do a:

    $ source conf/default.env
    $ python -m main.bin.useradd -u newuser -e john@gmail.com -n 'John Doe'
    *** creating user newuser
    $ python -m main.bin.postadd -e john@gmail.com -t 1 import/forum/rnseq-tool.txt
    *** adding 165, 104, john@gmail.com, RNA-SeqQC quality control

Versioning
-----------

BioStar uses the [common versioning nomenclature][versioning] that consists of
three numbers called `major.minor.revision` See the [NEWS.md][news] file for
release numbers and version migration information.

[versioning]: http://en.wikipedia.org/wiki/Software_versioning#Sequence-based_identifiers
[news]: https://github.com/ialbert/biostar-central/blob/master/NEWS.md

Testing
-------

API Testing can be initiated via the `biostar.sh` run manager:

    ./biostar.sh test

A `reports` directory will be created in the root directory that contains html
reports on the code coverage by the tests. View the `report/index.html` file.

Selenium tests can be run via:

    ./biostar.sh selenium
    
Please note that for this to work properly the python selenium library bindings
must be installed moreover the ``SELENIUM_TEST_LOGIN_TOKEN`` variable must be
set in your Django settings file. See the `conf/selenium.env` file:

    SELENIUM_TEST_LOGIN_TOKEN = "somepasswordgoeshere"

In addition both the test site and the command line above must make use of the
same settings file.

How the site works
-------------------

Posts may be formatted in [Markdown][markdown].

User reputation is a sum of all upvotes and accepted answers that a user
accumulates. Note that multiple answers may be accepted on a question, in
effect this provides the author of a question to reward twice the excellent
answers.

In Biostar there are four types of users: anonymous users, registered users,
moderators and administrators.

- *Anonymous users*: May browse all content of a site.
- *Registered users*: In addition to the privileges that anonymous users have
  registered users may create new posts if their reputation exceeds a limit
  (the default is zero), may vote and post answers and comments. 
- *Editors*: In addition to the privileges that registered users have Editors
  may edit, close and delete posts, edit user information (other than email)
  and may also suspend and reinstate users. All the actions of the Editors may be
  followed via the Moderator Log page (see About BioStar page for a link)
- *Administrators*: In addition to the privileges that moderators have
  administrators may promote/demote users from having moderator roles.
  Administrators also have access to the Django admin interface where they may
  perform more database actions than those offered via the BioStar interface..

Internationalization
--------------------

There are some dependencies that need to be installed (notable the Unix
`gettext` utility) to run the Django [makemessages][makemsg] command.

Template content needs to be tagged with the [Django Translation][trans]
framework.  The new message compilation then will be run via:

    ./biostar.sh messages

The settings file needs to specify the language (see the fileas named
`conf/ch.env` and `conf/ch.py`). For an example site in either Chinese run the
following:

    source conf/default.env
    source conf/ch.env
    ./biostar.sh run

To override what pages get loaded in add a new template directory in a location
and override the template loader order.  For an example see the `conf/ch.py`
settings file. In this example a different widget will be loaded from the
corresponding templates in `conf/custom-html/widgets/page.share.html`.

[makemsg]: https://docs.djangoproject.com/en/dev/ref/django-admin/#makemessages

[trans]: https://docs.djangoproject.com/en/dev/topics/i18n/translation/#internationalization-in-template-code


Content Persistence
-------------------

Content may be deleted (marked invisible to users) or destroyed (removed from
the database).

A post submitted for deletion will be destroyed only if the author requests the
deletion and the post does not have any follow-ups (answers/comments)
associated with it. Deleted top level posts are marked invisible to regular
users.

Code Layout
-----------

The Python code, templates, static content (css, images, javascript) and
default database are found in the *main* directory. There is partial datadump
of the existing BioStar content in the *import* folder. The *import* command
will load this data into the current database.

Other Libraries
---------------

Biostar is built with open source libraries. The following software packages
are used and if necessary included and distributed with BioStar:

* [Bootstrap][bootstrap] as a CSS framework
* [JQuery][jquery] for javascript programming
* [Less][less] used for syntactically awesome css
* [markitup][markitup] as rich text javascript editor. 
* [python-markdown2][markdown2] python library to convert [Markdown][markdown] to  HTML
* [docutils][docutils] is used to convert ReST_ to HTML
* [django_openid_auth][django_openid_auth] and [python_openid][python_openid]
    for openid authentication
* [whoosh][whoosh] provides fast full text searching
* [coverage][coverage] is used to measure code coverage during testing
* [prettify][prettify] is used for syntax highlighting
* [html5lib][html5lib] provides html parsing

[coverage]: http://pypi.python.org/pypi/coverage
[django_openid_auth]: https://launchpad.net/django-openid-auth
[python_openid]: http://pypi.python.org/pypi/python-openid/
[whoosh]: https://bitbucket.org/mchaput/whoosh/wiki/Home
[markdown2]: https://github.com/trentm/python-markdown2/
[django]: http://www.djangoproject.com/
[python]: http://www.python.org/
[jquery]: http://jquery.com/
[markitup]: http://markitup.jaysalvat.com/home/
[less]: http://lesscss.org/
[prettify]: http://code.google.com/p/google-code-prettify/
[bootstrap]: http://twitter.github.com/bootstrap/
[docutils]: http://docutils.sourceforge.net/docs/user/rst/quickstart.html
[rest]: http://docutils.sourceforge.net/docs/user/rst/quickstart.html
[markdown]: http://en.wikipedia.org/wiki/Markdown
[html5lib]: http://code.google.com/p/html5lib/

Data Migration from SE1
-----------------------

To load content from a StackExchange 1 XML datadump one needs to *migrate* the
data into the new schema. This is accomplished via the `migrate` command::

	$ ./biostar.sh migrate

This command in turn invokes the `main/migrate.py` script. Run this script
(note that the Django settings need to be properly set beforehand) with the -h
flag to see the flags it can take.::

    $ python -m main.migrate.py -h

.. note:: The `migrate` command used via the `biostar.sh` run manager makes use
  of an in memory database as specified in the `conf/memory.env` and
  `conf/memory.py` files.

The result of a data migration is a compressed json data fixture file that, in
turn, may be used via the *import* command::

    $ ./biostar.sh init import index

Account migration
-----------------

There is an automatic account migration based on the email provided by the
OpenID provider. Only the information from a subset of well known OpenID
providers are trusted enough to allow automatic account merging. Accepted
providers are: Google, Yahoo, Myopenid, LiveJournal, Blogspot, AOL, and
Wordpress. For other users manual migration of accounts will be required.
Users listed in the Django *ADMINS* settings will have full administration
privileges.

There is a postgresql database management script in `conf/pg-manager.sh` that
is used to facilitate data dumps and restoration.

Environment variables may be used to customize the behavior:

- `DJANGO_SETTINGS_MODULE`: the configuration module for Django
- `PYTHON`: the python executable that is to be invoked
- `FIXTURE`: output path to the (gzipped) file that will contain the data fixture
- `MIGRATE_PATH`: path to the directory that stores the StackExchange XML dump
- `MIGRATE_LIMIT`: the number of records to load from the XML dump

For a current Biostar run with about 4K users, 30K posts, 40K edits, 60K votes
generates about 300K database entries of various kinds. Data migration into a
fixture takes about 1 hour and 10Gb of RAM. This is an area that we could do a
lot better job (possibly orders of magnitude better).

The resulting data fixture is database independent and can now be loaded into
type database: sqlite, mysql, postgresql supported by Django. For example when
loading into postgresql it takes about 2 hours and 2Gb of RAM.

Note that the databases can be dumped and restored with far fewer resources.
Exporting directly into/from postgresql for example takes less than a few
minutes.

