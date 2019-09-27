## Documentation for Biostars Forum 

Here are steps to running and deploying the forum from scratch.


1. Create a virtual environment and clone most recent version of the forum.


2. Install dependencies. 


3. Run migrations and tests. 


4. Start a local server. 


### 1. Create a virtual environment and clone the repo.

Create a virtual environment by first downloading miniconda at https://docs.conda.io/en/latest/miniconda.html. 

After downloading the installation file, run the command : 

    $ bash installation_file.sh      

Once miniconda has been installed, create a virtual enviroment called `engine`.

    $ conda create -n engine python=3.7
    
Start the virtual enviorment by entering the command.

    $ conda activate engine
    
Clone or pull the most recent version of the forum. 

      (engine) $ git clone https://github.com/ialbert/biostar-central.git  # Clone a new branch
      
      (engine) $ git pull https://github.com/ialbert/biostar-central.git   # Pull into an exisiting 
      
      
### 2. Install dependencies. 

Activate the `engine` virtual enviorment.

    $ conda activate engine

Enter the `biostar-central` directory to install dependencies and requirements into the virtual enviorment.

Run `pip install -r conf/pip_requirements.txt ` to install all python requirements. 

Run `conda install -f conf/conda_requirements.txt` to install all anaconda requirements.

    (engine) user:~/biostar-central$ pip install -r conf/pip_requirements.txt      # Install python requirements.
    
    
    (engine) user:~/biostar-central$ conda install -f conf/conda_requirements.txt  # Install conda requirements.
    
 After dependencies have been installed, a migration needs to be made to create the database collect static files.
 
 ### 3. Run migrations and tests. 
 
Activate the `engine` virtual enviorment.

    $ conda activate engine
    
Run the command `make forum init` inside of the `biostar-central` directory to migrate and collect static files. 
 
    (engine) user:~/biostar-central$ make forum init    # Migrate and collect static files. 

A local database has now been created and the collected static files can be found in `biostar-central/export/static/`

Tests are run using : `make forum test`. 

    (engine) user:~/biostar-central$ make forum test    # Run tests. 
    
### 4. Start a local server 

Activate the `engine` virtual enviorment.

    $ conda activate engine
    
Enter the `biostar-central` directory and enter the command `make forum serve` to start a local server 

    (engine) user:~/biostar-central$ make forum serve    # Start local server

The site is available at http://127.0.0.1:8000/. 





 

