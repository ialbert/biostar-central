##Manage

There are a number of data management commands that come with Biostar.

### The biostar.sh manager

The **biostar.sh** shell command automatizes a number of commonly used tasks. Run it
with no parameters to get help on a typical usage::

    Usage:

      $ biostar.sh <command>

    Multiple commands may be used on the same line:

      $ biostar.sh init import run

    Commands:

      init      - initializes the database
      run       - runs the development server
      index     - initializes the search index
      test      - runs all tests
      env       - shows all customizable environment variables

      import    - imports the data fixture JSON_DATA_FIXTURE=import/default-fixture.json.gz
      dump      - dumps data as JSON_DATA_FIXTURE=import/default-fixture.json.gz
      delete    - removes the sqlite database DATABASE_NAME=biostar.db

      pg_drop        - drops postgres DATABASE_NAME=biostar.db
      pg_create      - creates postgres DATABASE_NAME=biostar.db
      pg_import f.gz - imports the gzipped filename into postgres DATABASE_NAME=biostar.db

    Use environment variables to customize settings. See the docs.

    DJANGO_SETTINGS_MODULE=biostar.settings.base

### Subcommands

In addition there are a  number of data management commands that are implemented for the each app.
Run::

    python manage.py help

And look for the output for the app ``[server]``, these commands will look like::

    [server]
        biostar_pg_dump
        delete_database
        import_biostar1
        import_mbox
        initialize_site
        prune_data
        usermod
        sqlfix
        sitemap
        user_crawl
        test_email
        test_task
        patch

You can run each of these subcommands with the `-h` flag to get more information on them.

### Command line tagging

There is a command line tool to perform content tagging based on a regular expression. The
invocation is::

	workon biostar
	source live/deploy.env
    python manage.py patch --tag "regexp:tag1,tag2,tag3"

Where the regular expression ``regexp`` will be searched against the content and when found matching
tags ``tag1``, ``tag2``, ``tag3`` will be applied. Example::

    python manage.py patch --tag "gff:gff,interval"

To detect what posts would be tagged but not actually perform the tagging pass the ``--dry`` command.
In that case only the post titles will be listed::

    python manage.py patch --tag "gff:gff,interval" --dry

This command will navigate through all questions in the database.

### Example commands

Frequently used commands::

    # Set the password for a user identified by their userid
    python manage.py usermod -u 2 -p abcde

    # Set the password for a user identified by their email
    python manage.py usermod -e foo@bar -p abcde

    # Rebuild the entire search index
    python manage.py rebuild_index

    # Reindex only what has changed in the last hour
    python manage.py update_index --age 1

    # Import 100 posts from a mbox file into biostar
    python manage.py import_mbox -f filename -l 100

    # Create a postgres database dump
    python manage.py biostar_pg_dump

### Merging Users

Create a space separated text file that contains the emails in the form::

    master_email alias_email1 alias_email2 ...

Then run the command::

	python manage.py patch --merge_users yourfile.txt

The command will move all content, votes and accounts associated with users identified by
the aliases into the master email. It then deletes the alias users. The effect of this
command cannot be reverted other than loading up a backup database dump.