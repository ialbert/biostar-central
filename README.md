# biostar-engine

## How to connect

    ssh www@planktontow.com
    
    ssh nba@planktontow.com
    
    ssh aswathy@planktontow.com
    
 
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
 
    make test_pull
    make main_pull
    make restart_nginx
    make restart_uwsgi
    
    
    