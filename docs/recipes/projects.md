
# Projects

The platform is project based. Each project is a collection of data, recipes and results.

Thus each project has three distinct sections:

1. Data (the input files)
2. Recipes (the code that processes the data)
3. Results (the directory that contains the resulting files of applying the recipe to data)

## Privacy

Within the management interface, all content is grouped into projects that may have public or private visibility. 
Content stored in public projects is readable without restrictions. 
Private projects will restrict access to members only.

1.  Public - viewable to everyone
2.  Private - viewable by collaborators
3.  Sharable - actively shared amongst a set of users

## Directory Structure 

Each project has a physical directory associated on the system located on the system. 


1. Projects directory
    - Each project has a directory with the data associated. 
2. Results directory
   - Location where the results of a recipe run are stored.
3. Table of contents directory
    - Contains table of content files for every data.

These directories all found in the media directory found in the `settings.py` under `MEDIA_ROOT`. The general structure is:

    media/
        projects/
           ...
        jobs/
           ...
        tocs/
            ... 
 

## Create a project

Creating a project can be done using the command line or web interface. 
Ensure you have a local server running or have access to a remote one when using  the web interface.

### Using command line

Use the management command `project` to create a project from command line.

    $ python manage.py project --help
    
    usage: manage.py project [-h] --pid PID [--name NAME] [--info INFO] [--public]
                         [--update] [--version] [-v {0,1,2,3}]
                         [--settings SETTINGS] [--pythonpath PYTHONPATH]
                         [--traceback] [--no-color] [--force-color]

    Creates a project.
    
    optional arguments:
      -h, --help            show this help message and exit
      --pid PID             Project id
      --name NAME           Project name
      --info INFO           File path or text of the project info
      --public              Makes project public
      --update              Updates the project selected by pid
      
      ...
   
      
To create a sample project, run the following:

    python manage.py project --name sample project --public --info "This is a sample" --uid sample
      

### Using web interface

Click on the `New Project` tab circled on the right. 

![](images/new-project.png)

This will bring you to a form to fill in the name, privacy, information, etc...


## Access

Before any actions a user takes on the platform, their access to the project is checked.

Access level are:

Read:

- Clone and copy recipe
- Read and copy data
- Read and copy results
- Create and edit their own recipes

Share:

- Includes all permission in `Read Access`
- Activated using a sharable project link


Write:

- Includes all permission in `Read Access`
- Upload new data 
- Delete objects
- Edit all recipes in projects
- Add or remove collaborators to the project 

Recipes can be misused so running them requires more privileges.

**Admins**,**staff**,and **trusted users** can run recipes with read or write access.


## Granting Access

Adding collaborators can be done using the command line or the interface. 
Ensure you have a local server running or have remote when using the web interface.


### Using command line

To add a user using command line use the managment command `add_user`:

    $ python manage.py add_user --help --fname user_file.csv

    usage: manage.py add_user [-h] [--fname FNAME] [--version] [-v {0,1,2,3}]
                              [--settings SETTINGS] [--pythonpath PYTHONPATH]
                              [--traceback] [--no-color] [--force-color]
    
    Add users
    
    optional arguments:
      -h, --help            show this help message and exit
      --fname FNAME         The CSV file with the users to be added. Must have
                            headers: Name, Email
    

With a sample csv file `user_list.csv` that looks like :

    user 1,  user1@email
    user 2,  user2@email


You can run the following command using the file:

    python manage.py add_user --fname user_list.csv
    
    
### Using web interface 

Click on the `Info` tab to view the project menu bar.

Click on the middle button labeled `Manage Access` 
![](images/manage-access-button.png)


 Search for users and select the access you would like. to give them
![](images/results.png)


# Data

Data may be uploaded or may be linked directly from a hard drive or from a mounted filesystem, thus avoiding copying and transferring large datasets over the web.  
For recipes that connect to the internet to download data, for example when downloading from the `Short Read Archive` the data does not need to be already present in the local server.

Notably the concept of “data” in our system is broader and more generic than on a typical file system. 
In our software “data” may be a single file, it may be a compressed archive containing several files or it may be a path to a directory that contains any number of files as well as other subdirectories. 
The programming interfaces for recipes can handle directories transparently and make it possible to run the same recipes that one would use for a single file on all files of an entire directory. 


## Data Types



## Upload data

Data can be added multiple ways.

Web interface options:
  - Upload a file
  - Write text
  - Import directory
  
Command line options:
  
  - Link a file directly from a hard drive
  
### Using command line

You can use the management command `data` to add or edit `Data` objects.

    $ python manage.py data --help 
    
    usage: manage.py data [-h] --pid PID [--did DID] [--update] [--path PATH]
                      [--text TEXT] [--name NAME] [--type TYPE] [--version]
                      [-v {0,1,2,3}] [--settings SETTINGS]
                      [--pythonpath PYTHONPATH] [--traceback] [--no-color]
                      [--force-color]

    Adds data to a project
    
    optional arguments:
      -h, --help            show this help message and exit
      --pid PID             Select project by unique uid
      --did DID             Select data by unique uid
      --update              Update the table of content for data --did.
      --path PATH           Path to the data to be added (file or directory)
      --text TEXT           A file containing the description of the data
      --name NAME           Sets the name of the data
      --type TYPE           Sets the type of the data
    
    

Link a sample directory, `/path/to/data/`q, to an existing project,  `project one`:

    
    $ python manage.py data --pid project one --path /path/to/data/ --name New data

### Using web interface


Open the `Data` tab inside of a project. 

Then click on the `New Data` tab on the right. 
![](images/new-data.png)

This opens another a form with two options.


1 . **Upload a file** -  Comes with size restrictions that come be found in the `settings.py` 

![](images/upload.png)


2 .  **Write text** - 10k character limit 
![](images/write.png)
    

## Import directory

Admin, staff, and trusted users can see an extra tab labeled `Import Data`















