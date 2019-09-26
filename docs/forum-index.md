## Documentation for Biostars Forum 


### Running and deploying the Biostar forum

Here are steps to running and deploying the forum from scratch.


1. Create a virtual environment and clone most recent version of the forum.


2. Install dependencies, run migrations and tests. 


3. Run a local server to make changes.


4. Deploy local changes to remote server.


#### 1. Create a virtual environment and clone the most recent version. 

Create a virtual environment by first download miniconda at https://docs.conda.io/en/latest/miniconda.html. 

After downloading the installation file, run the command : 

    $ bash installation_file.sh      

Once miniconda has been installed, create a virtual enviroment called `engine`.

    $ conda create -n engine python=3.7
    
Start the virtual enviorment by entering the command.

    $ conda activate engine
    
Clone or pull the most recent version of the forum. 

      $(engine) git clone https://github.com/ialbert/biostar-central.git  # Clone a new branch
      $(engine) git pull https://github.com/ialbert/biostar-central.git   # Pull into an exisiting 

#### 2. Clone or pull the most recent version

The first step is to pull the most recent version of the forum from github.
 
      $ git clone https://github.com/ialbert/biostar-central.git  # Clone a new branch
      $ git pull https://github.com/ialbert/biostar-central.git   # Pull into an exisiting 
      
 

