# The Biostar Engine

## Better solutions for an imperfect world.

![Biostar Engine Badge](biostar/engine/static/images/badge-engine.svg)

[python]: https://www.python.org/
[django]: https://www.djangoproject.com/

The Biostar Engine is a [Python 3.6][python] and [Django][django] based scientific data analysis oriented application server that can execute scripts over web while providing a graphical user interface for selecting the parameters of the scripts. The scripts may written in `bash`, may be a `Makefile`, `R` commands, basically any executable command.

We call the scripts as *recipes*.

The software has support for data storage and project management and this can be used as simple LIMS (Laboratory Information Management System)

## Installation

Our installation instructions rely on [conda][conda] though other alternatives would be also viable, virtual env, homebrew or even not using any environment altogether. 

We will demonstrate the installation using [conda][conda].

1\. Create a virtual environment

[conda]: https://conda.io/docs/

    conda create -y --name engine python=3.6
    source activate engine
    
2\. Clone the source server code and the recipe code:

    git clone git@github.com:biostars/biostar-engine.git
    git clone git@github.com:biostars/biostar-recipes.git
    
3.\ Install the python dependencies:

    # Switch to the engine directory.
    cd biostar-engine
    
    # Install server dependencies.
    pip install -r conf/python_requirements.txt
    
At this point the installation is complete.

4\. Start the server

All commands run through `make`. To initialize and run the test site use:

      make reset serve
   
Now visit <http://localhost:8000> to see your site running. 

The default admin email/password combination is: `1@lvh.me/testbuddy`.  You may change these in the settings.

## Bioinformatics environment

The following step is an optional, [bioconda][bioconda] specific requirement. 
Use it only if you also want to run all [bioconda][bioconda] tools that we have recipes for.
Make sure that you have set up [bioconda][bioconda] if you wish to run this!

    # Activate the environment.
    source activate engine
      
    # Switch to the engine directory.
    cd biostar-recipes
    
    # Install the conda dependencies.
    conda install --file conf/conda_requirements.txt

    # Add the recipes to the python path.
    python setup.py develop

[bioconda]: https://bioconda.github.io/

## Make commands

The Makefile included with the engine contains additional commands.

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


## Security considerations

**Note**: The site is designed to execute scripts on a remote server. In addition the site 
allows administrative users to change the content of these scripts. 

It is **extremely important** to restrict and guard access to all 
accounts with administrative privileges! 







    
