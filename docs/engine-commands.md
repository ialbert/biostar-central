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
      -h, --help  show this help message and exit
      -l, --load  Load to url from a directory. Load to database if --url is not
                  set.
      -d, --dump  Dump from a url to directory. Dump from database if --url is not
                  set.
      --pid PID   Project uid to load from or dump to.
      --url URL   Site url.
      --key KEY   API key. Required to access private projects.
      --dir DIR   Directory to store/load data from.
      
    # Dump project from remote url using --dump (-d) flag.
    $ python manage.py api project -d --pid tutorial --url https://www.bioinformatics.recipes

    # Load project to remote url using --load (-l) flag.
    $ python manage.py api project -l --pid tutorial --url https://www.bioinformatics.recipes --key=API_KEY

Data API:

    $ python manage.py api data --help
    optional arguments:
      -h, --help    show this help message and exit
      --path PATH   Path to data.
      --name NAME   Name of data.
      --pid PID     Project uid to create data in.
      --uid UID     Data uid to update.
      --update_toc  Update table of contents for data --uid.

    # Create new data object pointing to --path in project.
    $ python manage.py api data --pid tutorial --path DATA
    
    # Update the table of content for data --uid
    $ python manage.py api data --uid 32dd2w --update_toc
    
Recipe API:

    $ python manage.py api recipe --help
    optional arguments:
      -h, --help  show this help message and exit
      --uid UID   Recipe uid to load or dump.
      -l, --load  Load to url from a directory. Load to database if --url is not
                  set.
      -d, --dump  Dump from a url to directory. Dump from database if --url is not
                  set.
      --pid PID   Project uid to load from or dump to.
      --url URL   Site url.
      --key KEY   API key. Required to access private projects.
      --dir DIR   Directory to store/load data from.

    # Dump all recipes in project --pid from remote host.
    $ python manage.py api recipe -d --pid tutorial --url https://www.bioinformatics.recipes --key API_KEY
    
    # Load all recipes in project --pid to remote host.
    $ python manage.py api recipe -l --pid tutorial --url https://www.bioinformatics.recipes --key API_KEY
    
    
    
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
    