# Biostar Central


### Note on July 22nd, 2019
### We are migrating to an entirely new codebase!
### The current main site is run from branch: biostar2016

Software for better science.

**Biostar Central** is a [Python][python] and [Django][django] based collection of web applications that support science education.

Each web application may be deployed individually or in combination with all the others.

Applications:

- `forum` runs a Q&A forum inspired by StackOverflow see: [Biostars Q&A][biostars]
- `recipes` is a webapp that can run command line scripts via a web interface see: [Bioinformatics Recipes][recipes]

In addition `biostar-central` includes applications that are not meant to be used independently but in conjunction with other apps.

- `accounts` is a an application to manage user accounts (
- `emailer` is an application for sending notificationa and emails


[python]: https://www.python.org/
[django]: https://www.djangoproject.com/
[biostars]: https://www.biostars.org
[recipes]: https://www.bioinformatics.recipes
[handbook]: https://www.biostarhandbook.com
[conda]: https://conda.io/docs/

## Installation

The **Biostar Central**  requires [Python 3.6][python] or above.

Our installation instructions rely on [conda][conda] though other alternatives are equally viable.

Users may use `virtualenv`, `pipenv`, `homebrew`, `apt-get` etc, or they may opt to not use environment management tools at all.

    # Create a virtual environment.
    conda create -y --name engine python=3.6
    
    # Activate the python environment.
    source activate engine

    # Clone the source server code and the recipe code.
    git clone git@github.com:ialbert/biostar-central.git

    # Switch to the biostar-engine directory.
    cd biostar-central

    # Install server dependencies.
    pip install -r conf/pip_requirements.txt

The installation is now complete. All server management commands run through `make`. We provide different initializations depending on each app. For example to run a demonstration version of the forum app execute:

    make forum init serve

or to run the demonstration version of the  `recipe` app execute:

    make recipes init serve

In each case visit `localhost:8000` to view the site. To list all the additional commands that may be executed with


    python manage.py help --settings biostar.forum.settings

or for recipes:

    python manage.py help --settings biostar.recipes.settings

## Note

The software follows the recommended practices for developing [Django web applications][django] .

The [Django documentation][django] contains a wealth of information on the alternative ways to deploy the site on different infrasructure.


##  Additional confguration

To run bioinformatics oriented recipes further configuration of your environment may be necessary. You may perform these via the following instructions:

    conda config --add channels r
    conda config --add channels conda-forge
    conda config --add channels bioconda

    # Install the conda requirements.
    conda install --file conf/conda_requirements.txt

