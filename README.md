## Biostar Central

### Software for better science.

**Biostar Central** is a [Python][python] and [Django][django] based collection of web applications that support scientific practice and education.

The goal of the project is to produce software with straightforward installation and minimal dependencies that works on any computing platform that supports Python. For each app the philosophy is that of decentralization and self hosting. We write code to allow others to recreate the same services that we run.

Each web application may be deployed individually or in combination with  the others. The following applications are currently feature complete:

- `recipes` a web app that can runs data analysis scripts via a web interface, see: [Bioinformatics Recipes][recipes]
- `forum` a web app that runs a Q&A forum inspired by StackOverflow, see: [Biostars Q&A][biostars]

Note: The public [biostars.org][biostars] site runs the code from the `biostar2016` branch. The  master branch is the development version that we will migrate the main server to in the future.


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
    pip install -r conf/pip_requirements.txt

The installation is now complete.

All server management commands are run through `make`.

## Testing the code

To test the `recipes` app run:

    make recipes test

To test the `forum` app run:

    make forum test

## Running the demo server

To run the demonstration version of the `recipes` app execute:

    make recipes_demo

To run a demonstration version of the `forum` app execute:

    make forum_demo

Visit <http://127.0.0.1:8000/> to view the site.

## Default users

All users listed in the `ADMINS` attribute of the Django `settings.py` module when the site is initialized. The
`DEFAULT_ADMIN_PASSWORD` attribute will be set as their password. By default the value for both is:

    admin@localhost

Use the username and password above to log into the site as an administrator.

## Executing multiple tasks

The Makefile has several tasks that demonstrate the commands that may be chained together:

Initialize and run an empty `recipes` app.

    make recipes init serve

Initialize and run an empty `forum` app.

    make forum init serve

## Documentation

Additional documentation for:

* [Documentation for the `recipes` app](https://bioinformatics-recipes.readthedocs.io/en/latest/index.html)
* [Documentation for the `forum` app](docs/forum/forum-index.md)


