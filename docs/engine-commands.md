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


Import data from a site using API:
    
    python manage.py api_import
    
    arguments:
        --key = API key only required for private projects.
        --base = Base site url to do a reverse look up of api urls. Defaults to url from settings.
        --output = Directory to download data to. Defaults to API_DUMP in settings.py.

**Import data: Example Usage**
        
    # Download local data into export/media/test
    python manage.py api_import --output=export/media/test
    
    # Download remote data 
    python manage.py api_import --base="https://www.bioinformatics.recipes/" --key=API_KEY
    
    
Export data to a site using API:


     python manage.py api_export --key=API_KEY
     
     arguments:
        --key = Required API key.
        --base = Base site url to do a reverse look up of api urls. Defaults to url from settings.
        --data = Base data directory to export multiple project data from. 
                 Its subdirectories are expected to be: /project/recipe.
        --project = Full path to single project directory to crawl and export recipes from.
        --recipe = Full path to single recipe dir to export from.
    

**Export data: Example Usage**
    
    # Upload multiple projects from a base dir
    python manage.py api_export --key=API_KEY --base=base/project1/recipe2 
    
    # Upload multiple recipes from a project dir
    python manage.py api_export --key=API_KEY --project=project3/
    
    # Upload a single recipe
    python manage.py api_export --key=API_KEY --recipe=recipe4/
 
    
    