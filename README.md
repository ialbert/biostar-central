## Biostar Central


---

### Note on July 22nd, 2019

We are migrating to an entirely new codebase!

The current main site is run from the `biostar2016` branch. We plan to switch over to the new codebase in the near future.

---

#### Software for better science.

**Biostar Central** is a [Python][python] and [Django][django] based collection of web applications that support science education.

Each web application may be deployed individually or in combination with all the others.

Application that currently operational:

- `recipes` is an app that runs scripts via a web interface. See: [Bioinformatics Recipes][recipes]
- `forum` is an app that runs a Q&A forum inspired by StackOverflow. See: [Biostars Q&A][biostars]

In addition `biostar-central` includes applications that are not meant to be used independently but in conjunction with other apps.

- `accounts` is an app that manages user accounts
- `emailer` is an app that can send emails


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

The installation is now complete. To test the code run:

    make forum test
    make recipe test

## Running the code

All server management commands run through `make`. We provide different initializations depending on each app. For example to run a demonstration version of the `forum` app execute:

    make forum init serve

To run the demonstration version of the `recipes` app execute:

    make recipes init serve

In each case visit `http://localhost:8000` to view the site.

The admin username and password are set via the `ADMINS` and the `ADMIN_PASSWORD` settings in `biostar/accounts/settings.py`. By default both the admin login name and the default admin password are set to

    admin@localhost

**Note**: These settings must be changed on a publicly facing site!

## Documentation

Additional documentation in [docs/index.md](docs/index.md)

## Django configuration

The software follows the recommended practices for developing [Django web applications][django] .

The [Django documentation][django] contains a wealth of information on the alternative ways to deploy the site on different infrastructure.


##  Infrastructure configuration

To run bioinformatics oriented software via the recipes additional configuration of your environment may be necessary. For example we use:

    conda config --add channels r
    conda config --add channels conda-forge
    conda config --add channels bioconda

    # Install the conda requirements.
    conda install --file conf/conda_requirements.txt

