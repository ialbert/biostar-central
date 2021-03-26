## Biostar Central

### Software for better science.

**Biostar Central** is a [Python][python] and [Django][django] based collection of web applications that support scientific practice and education.

The goal of the project is to produce software with straightforward installation and minimal dependencies that works on any computing platform that supports Python. For each app the philosophy is that of decentralization and self hosting. We write our code to allow others to recreate the same services that we run.

Each web application may be deployed individually or in combination with  the others. The following applications are currently feature complete:

- `forum` a web app that runs a Q&A forum inspired by StackOverflow, see: [Biostars Q&A][biostars]
- `recipes` a web app that can runs data analysis scripts, see: [Bioinformatics Recipes][recipes]

[python]: https://www.python.org/
[django]: https://www.djangoproject.com/
[biostars]: https://www.biostars.org
[recipes]: https://www.bioinformatics.recipes
[handbook]: https://www.biostarhandbook.com
[conda]: https://conda.io/docs/

## Installation

The code in **Biostar Central**  requires [Python 3.6][python] or above.

Our installation instructions rely on [conda][conda] though other alternatives for managing python environments are equally viable.

    # Create a virtual environment.
    conda create -y --name engine python=3.6
    
    # Activate the python environment.
    conda activate engine

    # Clone the source server code and the recipe code.
    git clone https://github.com/ialbert/biostar-central.git

    # Switch to the biostar-engine directory.
    cd biostar-central

    # Install server dependencies.
    pip install -r conf/requirements.txt

The installation is now complete.

All server management commands are run through `make`.

## Demo server

To run the demonstration version of the `recipes` app execute:

    make recipes demo

To run a demonstration version of the `forum` app execute:

    make forum demo

Visit <http://127.0.0.1:8000/> to view the site.

## Default users

All users listed in the `ADMINS` attribute of the Django `settings.py` 
module will gain administritave privileges when the site is initialized. The
`DEFAULT_ADMIN_PASSWORD` attribute will be set as the default admin password. 
By default the value for both is:

    admin@localhost

Use this username and password combination to log into the site as an administrator. Change the `DEFAULT_ADMIN_PASSWORD` for public facing installations.

## Running the site

The Makefile has several tasks that demonstrate the commands that may be run. Typically a series of `make` tasks may be run. For example

Initialize and run a new `recipes` app.

    make recipes serve

Initialize and run a demo version of the recipes app:

    make recipes demo
    
Initialize and run an new `forum` app.

    make forum serve

## Valid tasks

- `forum`, `recipes`: selects the app to run. It must be the first task in the list.
- `demo`: runs a demonstration version
- `init`: initializes the database schema
- `serve`: runs the app on the default port (`localhost:8080`)
- `save`: saves the current database content as a JSON fixture file
- `load`: loads the last database save file from a JSON fixture
- `reset`: resets the database (deletes all database content)
- `hard_reset`: runs `reset` then deletes all files in the media/spool folder

## Testing

To run all tests type:

    make test
    
To test the `recipes` app run:

    make recipes test

To test the `forum` app run:

    make forum test

## Documentation

Additional documentation for:

* [Documentation for the `recipes` app](https://bioinformatics-recipes.readthedocs.io/en/latest/index.html)
* [Documentation for the `forum` app](docs/forum/forum-index.md)


