# Scripts available to the forum app


## spam.sh

Clear spam from spam index

    # Clear spam from search index.
    python manage.py index --clear_spam


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


    