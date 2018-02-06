# The Biostar Engine

### Better solutions for an imperfect world.

The Biostar Engine is [Python 3.6][python] and [Django][django]based scientific data analysis oriented application server.

[python]: https://www.python.org/
[django]: https://www.djangoproject.com/

![Biostar Engine Badge](biostar/engine/static/images/badge-engine.svg)

The Biostar Engine allows executing scripts over the web via graphical user interfaces.
The scripts may written in `bash`, `R` or any other a scripting enabled language.

We call these scripts as "recipes".

In addition the software has support for data storage and project management. 
The Biostar Engine can be used as simple LIMS (Laboratory Information Management System)

## Installation

We recommend using [conda][conda]. For simplicity our installation 
instructions rely on [conda][conda] though other alternatives would be also viable, virtual env,
homebrew etc. Create a virtual environment both on your system and on the remote site:

[conda]: https://conda.io/docs/

    conda create -y --name engine python=3.6
    source activate engine
    
Clone the source server code and the recipe code:

    git clone git@github.com:biostars/biostar-engine.git
    git clone git@github.com:biostars/biostar-recipes.git
    
Install python dependencies:

    # Switch to the engine directory.
    cd biostar-engine
    
    # Install server dependencies.
    pip install -r conf/python_requirements.txt
    
    # The current package has packages that are needed at the command line.
    python setup.py develop
    
The following step is an optional, [bioconda][bioconda] specific requirement. 
Use it only if you also want to run all [bioconda][bioconda] tools that we have recipes for.
Make sure that you have set up [bioconda][bioconda] if you wish to run this!

[bioconda]: https://bioconda.github.io/

    conda install --file conf/conda_requirements.txt
    
## Quick start

All commands run through `make`. To initialize and run the test site use:

    make reset serve
    
Now visit <http://localhost:8000> to see your site running.

The default admin email/password combination is: `1@lvh.me/testbuddy`. 
You may change these in the settings.

## Valid commands

Re-initialize the database:

    make reset 
 
Serve the current site:

    make serve

Get the data that are used in the demonstration recipes.

    make data
            
Loads example recipes from the `biostar-recipe` repository that you cloned in the setup.

    make recipes

Run all tests:

    make test
        
## Deployment

The site is built with Django hence the official Django documentation applies.

* <https://docs.djangoproject.com/>

The software also supports `uwsgi` as the runtime architecture. When deplying through 
`uwsgi` jobs are queued and run automatically through the `uwsgi` spooler. See the `uwsgi` documentation 
for details on how to control that process.

* <https://uwsgi-docs.readthedocs.io/en/latest/>

The jobs may also be started as commands. See the `job` command for details:

    python manage.py job --help
    
For example

    python manage.py job --next
    
will execute the first job in the queue.

## Recipes

Recipes are stored and distributed from a separate repository at:

* <https://github.com/biostars/biostar-recipes>


## Security warning!

**Note**: The site is designed to execute scripts on a remote server. In addition it 
allows administrative users to change the content of these scripts. 
It is **extremely important** to properly restrict and guard access to
accounts with administrative privileges! 







    
