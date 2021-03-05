# Scripts available to the forum app


## celery_beat.sh

Start celery beat 

    $CELERY -A $APP beat -l info -f $LOGFILE


## celery_worker.sh

Start celery worker

    $CELERY -A $APP worker -l info --maxtasksperchild $MAX_TASK --concurrency $NUM_WORKERS -f $LOGFILE


## index.sh

Start search indexing based on batch size

    # Add 5000 posts to search index every 3 minutes
    python manage.py index --index ${BATCH_SIZE} --report
    

## planet.sh

Update planet app from

    # Update latest five entries for each planet blog.
    python manage.py planet --update ${UPDATE_COUNT}


## transfer.sh

Transfer from old data base to new ( several setting )

    
    #  Migrate newly created biostar next database
    python manage.py migrate --settings ${TRANSFER_SETTINGS_MODULE}


    # Transfer the old database into new database.
    python manage.py transfer --limit $LIMIT --settings ${TRANSFER_SETTINGS_MODULE}


    