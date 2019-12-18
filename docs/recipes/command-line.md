# Commands

Site administrator with shell access to the server can use these commands to interact with the recipes platform in an automated fashion.



## Create a Project

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

Note: The owner of any project created from command line is an first admin user.

      
To create a sample project, run the command:

    python manage.py project --name sample project --public --info "This is a sample" --pid sample
 

## Granting Access

Adding collaborators can be done using the command line or the interface. 

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
    
    
## Upload Data

Command line options:
  
  - Link a file directly from a hard drive
  
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
    
    

Link a sample directory, `/path/to/data/`, to an existing project with the uid  `project_one`:

    
    $ python manage.py data --pid project_one --path /path/to/data/ --name New data

## Create a Recipe 

Creating a recipe can be by directly upload json and script template to a given recipe.

Use the `recipe` management command to directly add to a project.

    $ python manage.py recipe --help
    
    usage: manage.py recipe [-h] --pid PID --rid RID [--json JSON]
                        [--template TEMPLATE] [--info INFO] [--name NAME]
                        [--image IMAGE] [--update] [--version] [-v {0,1,2,3}]
                        [--settings SETTINGS] [--pythonpath PYTHONPATH]
                        [--traceback] [--no-color] [--force-color]

    Adds recipe to a project

    optional arguments:
      -h, --help            show this help message and exit
      --pid PID             Project id.
      --rid RID             Recipe id.
      --json JSON           Recipe json path.
      --template TEMPLATE   Recipe template path (optional)
      --info INFO           Recipe description (optional)
      --name NAME           Recipe name
      --image IMAGE         Recipe image path
      --update              Updates the recipe


For example, the command below would add a recipe named `New recipe` to project with uid `1`.

    python manage.py recipe --pid 1 --name New recipe --json < interface file > --template < script template > 

