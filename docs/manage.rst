Manage
======

There are a number of data management commands that come with Biostar.

The biostar.sh manager
----------------------

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

Subcommands
-----------

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

You can run each of these subcommands with the `-h` flag to get more information on them.

Example commands
----------------

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