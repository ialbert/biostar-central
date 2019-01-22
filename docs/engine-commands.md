## Additional commands

The Makefile included with the engine contains additional commands.

Test the software:

    make test

Re-initialize the database:

    make reset

Serve the current site:

    make serve

Initialize the example recipes from the `biostar-recipe` repository.

    make recipes

Run all tests:

    make test

Back up the data

    make backup


Dump and load recipes from/to remote site:
    
    python manage.py api 
    
      -l, --load            Load data to url from a directory. Load recipes to
                            database if --url is not set.
      -d, --dump            Dump recipes from a url to directory. Dump recipes
                            from database if --url is not set.
      --url URL             Site url.
      --key KEY             API key required to get private projects.
      --dir DIR             Directory to store/load data from.
      --pid PID             Project uid to load or dump.
      --rid RID             Recipe uid to load or dump.
      
**Example Usage**
    
    # Dump recipe from remote site 
    python manage.py api -d --pid=cookbook --url="https://www.bioinformatics.recipes/" 
    
    # Load recipe into 
    python manage.py api_import --base="https://www.bioinformatics.recipes/" --key=API_KEY
    
    
    