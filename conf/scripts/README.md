#Scripts used in production servers for certain tasks


## db_backup.sh

Run pg_dump tasks behind to dump the most recent copy of database.


    # pg_dump the database
    python manage.py tasks --action pg_dump --outdir export/backup


## cleanup.sh

Run clean up management task

    python manage.py cleanup

## migrate.sh

Run collect static and migration

    # Migrate the server.
    python manage.py migrate
    
    # Collect static files
    python manage.py collectstatic --noinput

## reset.sh

Reset the database (if sql), cached static files, and logs.

    #python manage.py flush --noinput
    make reset