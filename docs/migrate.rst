Migrating from Biostar 1.X
==========================

Due to the complete rework there is no database schema migration.

Instead users of
Biostar 1 site are expected to export their data with a script provided in Biostar 1
then import it with a management command provided with Biostar 2.

The migration will take the following steps:

1. Set the ``BIOSTAR_MIGRATE_DIR`` environment variable to point to a work directory that
   will hold the temporary data, for example  ``export BIOSTAR_MIGRATE_DIR="~/tmp/biostar_export"``

2. Load the environment variables for the Biostar 1 site
   then run ``python -m main.bin.export -u -p -v``. This will dump the contents of the site
   into the directory that ``BIOSTAR_MIGRATE_DIR`` points to.

3. Load the environment variables for you Biostar 2 site then run the
   ``./biostar.sh import_biostar1`` command.

Some caveats, depending how you set the variables you may need to be located in
the root of your site. This applies for the default settings that both sites come
with, as the root is determined relative to the directory that the command is run in.