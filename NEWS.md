Tagged Production Releases
==========================

BioStar 1.2.3 on Sept 5, 2012: added shortcut commands

BioStar 1.2.2 on Aug 24, 2012: changed database schema to support bookmark
    counts in search queries, enabled sticky posts

BioStar 1.1.0 on Aug 12, 2012: changed site layout, database schema supports
    post view counting

Biostar 1.0.0 on April 7, 2012: migrated from the StackExchange 1 datadumps
    into to the www.biostars.org domain

Versioning
==========

BioStar uses the [common versioning nomenclature][versioning] that consists of
three numbers called `major.minor.revision`

Versions will be presented as three numbers called major.minor.revision.
Database dumps can only be loaded into the same major + minor version that they
were created under.

The version number is located in the file `main/server/__init__.py` can be
automatically extracted via the command:

    python -c "import main.server; print main.server.VERSION"
    
Changes in either major or minor revisions will require a database migration.
This can be accomplished by invoking:

    python manage.py migrate main.server 

New installations of BioStar are automatically migrated to the latest version. 

[versioning]: http://en.wikipedia.org/wiki/Software_versioning#Sequence-based_identifiers

