# biostar-engine

## How to connect

    ssh www@metabarcode.com
    
    ssh nba@metabarcode.com
    
    ssh aswathy@metabarcode.com
    
    
## Server location

    /export/www
    
## Install

Create a virtual environment both on your system and on the remote site:

    conda create -y --name engine python=3.6
    source activate engine
    
Switch to the server installatio directory to install conda and python dependecies:

    # Skip this step if you don't wish to run the tools that come with the engine.
    conda install --file conf/conda_requirements.txt
    
    # Install server dependencies.
    pip install -r conf/python_requirements.txt
    
    # The current package has packages that are needed at the command line.
    python setup.py develop
    

## Usage

import data into local machine

    make testdata

Creates a brand new database

    make reset 

Migrate and initialize the data:

    make init
   
Server local development:
   
    make serve

## Fabric scripts

Fabric works with python 2 only.
     
    conda create -y --name py2 python=2
    source activate py2
    pip install fabric
    
Make rules:

    make deploy_test
    make deploy_main
    make restart_nginx
    
    # Needs sudo privi for this command thats why www does not work
    make USER=nba restart_django


    
    
