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


Project API:

    $ python manage.py api project --help
    
    optional arguments:
      -h, --help            show this help message and exit
      -l, --load            Load to url from a directory. Load to database if
                            --url is not set.
      -d, --dump            Dump from a url to directory. Dump from database if
                            --url is not set.
      -u, --url_from_json   Extract url from conf file instead of --url. Only
                            works when --load(-l) flag is set.
      --url URL             Site url.
      --key KEY             API key. Required to access private projects.
      -a, --add_data        Add data found in conf file to local database.
      --privacy PRIVACY     Privacy of project, only used when creating.
      --uid UID             Project uid to load from or dump to.
      --dir DIR             Directory to store/load project from.
      --list                Show a project list.
      --data_root DATA_ROOT
                            Root directory to data found in conf file.
      --json JSON           JSON file path relative to --dir to get conf from ONLY
                            when --load flag is set.


    # Dump project from remote url using --dump (-d) flag.
    $ python manage.py api project -d --uid tutorial --url https://www.bioinformatics.recipes

    # Load project to remote url in json file using --load (-l)  
    # and --url_from_json (-u) flag.
    
    $ python manage.py api project -l -u --json tutorial.hjson --key API_KEY
    
Data API:

    $ python manage.py api data --help
    optional arguments:
      -h, --help    show this help message and exit
      --path PATH   Path to data.
      --pid PID     Project uid to create data in.
      --uid UID     Data uid to load or update.
      --text TEXT   A file containing the description of the data
      --name NAME   Sets the name of the data
      --type TYPE   Sets the type of the data
      --list        Show a data list.
      --update_toc  Update table of contents for data --uid.
      
    # Create new data object pointing to --path in project.
    $ python manage.py api data --pid tutorial --path DATA
    
    # Update the table of content for data --uid
    $ python manage.py api data --uid 32dd2w --update_toc
    
Recipe API:

    $ python manage.py api recipe --help
    optional arguments:
      -h, --help           show this help message and exit
      -l, --load           Load to url from a directory. Load to database if --url
                           is not set.
      -d, --dump           Dump from a url to directory. Dump from database if
                           --url is not set.
      -u, --url_from_json  Extract url from conf file instead of --url. Only works
                           when --load(-l) flag is set.
      --url URL            Site url.
      --key KEY            API key. Required to access private projects.
      --jobs               Also creates a queued job for the recipe
      --uid UID            Recipe uid to load or dump.
      --pid PID            Project uid to load from or dump to.
      --dir DIR            Directory to store/load recipe from.
      --list               Show a recipe list.
      --json JSON          JSON file path relative to --dir to get conf from ONLY
                           when --load flag is set.

    # Dump all recipes in project --pid from remote host.
    $ python manage.py api recipe -d --pid tutorial --dir tutorial --url https://www.bioinformatics.recipes 
    
    # Load all recipes in project directory --dir to remote host.
    $ python manage.py api recipe -l -u --dir tutorial --key API_KEY
    
    
    
Job API:

    $ python manage.py api job --help
    
    optional arguments:
      -h, --help            show this help message and exit
      --next                Runs the oldest queued job
      --id ID               Runs job specified by id.
      --uid UID             Runs job specified by uid.
      --show_script         Shows the script.
      --show_json           Shows the JSON for the job.
      --show_template       Shows the template for the job.
      --show_command        Shows the command executed for the job.
      --use_json USE_JSON   Override the JSON with this file.
      --use_template USE_TEMPLATE
                            Override the TEMPLATE with this file.
      --list                Show a job list
      
      
    # Run a job overriding json
    $ python manage.py api job --uid 332eqwd --use_json JSON_FILE
    