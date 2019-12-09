# General Concepts


## What is a recipe?

Each recipe is built from two ingredients:

1. The interface specification file.
2. The template specification file.

The **interface** will specify the value of the parameters that get substituted into the **template**.

The **template** contains the commands that need to be executed. The **template** will have
placeholders for the parameter values that the user will need to enter in the interface.

The interface + template will generate a script that the site can execute.

The software will generate an web interface for each parameter specified in the interface. It is this interface where users are able to select the values that their recipe needs to operate.


A recipe consists of a "JSON definition file" and a "script template".

The simplest JSON definition file is

    {}

A simple script template can be completely empty or might contain just:

    echo 'Hello World!'
    
 
## Minimum installation requirements for web application


Bioinformatics Recipes is a [Python](<http://www.python.org/>) and
[Django](<http://www.djangoproject.com/>) based data analysis software licensed under the *MIT Open Source License*.
It is a simple, generic, flexible and extensible analytic framework.

Requirements:
- **Python 3.7** or greater
- **Anaconda**
- **Django 2.6** or greater

The sourcecode can be obtained via::

    git clone https://github.com/ialbert/biostar-central.git
    
1 . Activate the virtual environment.

 
    $ conda activate < virtual environment >

2 . Enter the `biostar-central` directory to install dependencies and requirements into the virtual environment.
    
    
3 . Execute the following to install python requirements: 


    $ pip install -r conf/pip_requirements.txt      # Install python requirements.
    
4 . Execute the following to install all anaconda requirements:
    
    conda install --file conf/conda_requirements.txt  # Install conda requirements.
    
## Run a local server

A local server can be user 

## Who can execute recipes?

Only **admins**,**staff**,and **trusted users** can execute recipes. 

Admins, staff, and trusted users also need read or write access to the recipe. 

The different access levels are `Share Access`, `Read Access`, `Write Access`.


Read/Share Access:

- Clone/copy recipe
- Read/copy data
- Read/copy results
- Create/edit their own recipes

Owner/Write Access:

- Clone/copy/edit/create recipes
- Read/copy/upload data
- Read/copy/ results
- Edit all recipes in projects
- Add collaborators to the project 

       
## How do I run recipes?

The local web server on `localhost:8000` you can log in using `admin@localhost` 


## Creating new recipes

## Where do recipes actually run when executed?


## How do we standardize the way recipes amongst multiple users?


**Cloning recipes**

## What conventions should be followed when creating recipes


## Documenting recipes


## What is the amount of provenance (historical record) that recipes produce?


## Presenting the recipe provenance and result.


## How does Bioinformatics Recipes differ from Galaxy


### Advantages 


###  Functional differences seen by users


### Performance comparison 


## Scaling to larger users




