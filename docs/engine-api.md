# API Methods and Commands

## Methods

### Project list 
    GET /api/project/list/

List of projects in a tab delimited fashion with three columns: **ID** , **Name**, and **Privacy**

##### Example
[/api/project/list/](https://www.bioinformatics.recipes/api/project/list/)

    tutorial	Recipe Tutorials	Public
    cookbook	Bioinformatics Cookbook	Public

###Project Information

    GET /api/project/info/{id}/
    PUT /api/project/info/{id}/
   
##### Parameters
* _id_: unique project ID

##### Example
[/api/project/tutorials/](https://www.bioinformatics.recipes/api/project/tutorials/)

    {
        settings:
        {
            uid: tutorial
            name: Recipe Tutorials
            image: tutorial.png
            privacy: Public
            help:
                '''
                This project contains simple analyses that demonstrate the process
                of creating a **recipe**.
    
                Follow the **instructions,** investigate the **data** and **recipe code**
                to gain a deeper understanding of how recipes work.
    
                Read the step by step instructions in the [How to write recipes](https://github.com/biostars/biostar-recipes/blob/master/docs/how-to-write-recipes.md).
                '''
            id: 2
            project_uid: tutorial
            url: http://localhost:8000
        }
        recipes:
        [
            empty
            makefile
            starter
            interface
            hello-world
            environment
            rscript
        ]
    }

### Project Image

    GET /api/project/image/{id}/
    PUT /api/project/image/{id}/
  

##### Example

[/api/project/image/tutorials/](tutorial.png)



##### Parameters
* _id_: unique project ID

### Recipe list
    GET /api/recipe/list/{id}/
    
List of recipes for project in a tab delimited fashion with two columns: **ID** and **Name**

##### Parameters
* _id_: unique project ID

##### Example
[/api/recipe/list/tutorial/](https://www.bioinformatics.recipes/api/recipe/list/tutorial/)

    interface	Interface Elements
    makefile	Makefile Example
    environment	Environment Check

### Recipe Json

    GET /api/recipe/json/{id}/
    PUT /api/recipe/json/{id}/

Recipe JSON used to generate interface

##### Parameters
* _id_: Unique recipe ID

##### Fields in response 
Fields associated with the recipe JSON

##### Example
[/api/recipe/json/93412cee/](https://www.bioinformatics.recipes/api/recipe/93412cee/json/)

    {

      size: {
      
        label: Sample size
        display: INTEGER
        range: [ 1, 10000]
        value: 1000
      }

      settings: {
        name: R Script
        summary: This recipe demonstrates the use of an R script as a recipe.
        image: rscript.jpg
         execute: {
                filename: "recipe.r"
                command: "Rscript --vanilla recipe.r"
            }
        help:
        '''
        # Help

        This recipe is a demonstration of using an R script as a recipe.

        It creates a barplot of the `sample` function.

        '''
      }

    }

### Recipe Template

    GET /api/recipe/template/{id}/
    PUT /api/recipe/template/{id}/

Recipe template executed during analysis.

#### Parameters
* _id_: Unique recipe ID

#### Example
[/api/recipe/template/93412cee/](https://www.bioinformatics.recipes/api/recipe/93412cee/template/)

    # Set graphics device to PNG.
    fname = 'plot.png'
    png(fname)

    # Generate the sample.
    data = sample(1:3, size={{size.value}}, replace=TRUE, prob=c(.30,.60,.10))

    # Turn it into a table.
    data = table(data)

    # Generate the barplot.
    barplot(data)

    # Tell the user what happened.
    sprintf("Saved plot into file: %s", fname)


##Commands 

Recipe API:

    $ python manage.py api recipe --help
    optional arguments:
      -h, --help           show this help message and exit
      -l, --load           Load to url from a directory. Load to database if --url
                           is not set.
      -d, --dump           Dump from a url to directory. Dump from database if
                           --url is not set.
      -u, --url_from_json  Extract url from conf file instead of --url. Only works
                           when --load(-l) flag is set.
      --url URL            Site url.
      --key KEY            API key. Required to access private projects.
      --jobs               Also creates a queued job for the recipe
      --uid UID            Recipe uid to load or dump.
      --pid PID            Project uid to load from or dump to.
      --dir DIR            Directory to store/load recipe from.
      --list               Show a recipe list.
      --json JSON          JSON file path relative to --dir to get conf from ONLY
                           when --load flag is set.

    # Dump all recipes in project --pid from remote host.
    $ python manage.py api recipe -d --pid tutorial --dir tutorial --url https://www.bioinformatics.recipes 
    
    # Load all recipes in project directory --dir to remote host.
    $ python manage.py api recipe -l -u --dir tutorial --key API_KEY
    
    
    
Project API:

    $ python manage.py api project --help
    
    optional arguments:
      -h, --help            show this help message and exit
      -l, --load            Load to url from a directory. Load to database if
                            --url is not set.
      -d, --dump            Dump from a url to directory. Dump from database if
                            --url is not set.
      -u, --url_from_json   Extract url from conf file instead of --url. Only
                            works when --load(-l) flag is set.
      --url URL             Site url.
      --key KEY             API key. Required to access private projects.
      -a, --add_data        Add data found in conf file to local database.
      --privacy PRIVACY     Privacy of project, only used when creating.
      --uid UID             Project uid to load from or dump to.
      --dir DIR             Directory to store/load project from.
      --list                Show a project list.
      --data_root DATA_ROOT
                            Root directory to data found in conf file.
      --json JSON           JSON file path relative to --dir to get conf from ONLY
                            when --load flag is set.


    # Dump project from remote url using --dump (-d) flag.
    $ python manage.py api project -d --uid tutorial --url https://www.bioinformatics.recipes

    # Load project to remote url in json file using --load (-l)  
    # and --url_from_json (-u) flag.
    
    $ python manage.py api project -l -u --json tutorial.hjson --key API_KEY
    