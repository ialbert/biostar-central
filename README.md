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
    pip install -r requirements.txt

## Usage

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

    # example
    CURRENT=nba

    make deploy_test
    make deploy_main
    make restart_nginx
    
    # Needs sudo privi for this command thats why www does not work
    make USER=$CURRENT restart_django


    
    
