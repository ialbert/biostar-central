# API Methods and Commands



## Commands 

Dump API:

    $ python manage.py api dump --help
    
    optional arguments:
      -h, --help     show this help message and exit
      -r, --recipes  Load recipes of --pid
      --url URL      Site url.
      --key KEY      API key. Required to access private projects.
      --rid RID      Recipe uid to dump.
      --pid PID      Project uid to load dump.
      --dir DIR      Directory to store from.

    # Dump project from remote url (--url).
    $ python manage.py api dump --pid tutorial --dir tmp/remote/projects/ --url URL

    # Dump recipes from remote url ( --url) 
    $ python manage.py api dump --pid tutorial --dir tmp/remote/recipes/ --url URL --recipes
    
    
Load API:
    
    $ python manage.py api load --help
    optional arguments:
      -h, --help            show this help message and exit
      -u, --url_from_json   Extract url from conf file instead of --url.
      -r, --recipes         Load recipes of --pid
      -d, --data            Load data of --pid to local database.
      --url URL             Site url.
      --key KEY             API key. Required to access private projects.
      --add_data            Add data found in --json to --pid.
      --jobs                Also creates a queued job for the recipe
      --rid RID             Recipe uid to load.
      --pid PID             Project uid to load from or dump to.
      --did DID             Data uid to load or update.
      --dir DIR             Base directory to store/load from.
      --path PATH           Path to data.
      --text TEXT           A file containing the description of the data
      --name NAME           Sets the name of the data
      --type TYPE           Sets the type of the data
      --update_toc          Update table of contents for data --uid.
      --data_root DATA_ROOT
                            Root directory to data found in conf file when loading
                            project.
      --json JSON           JSON file path relative to --dir.
                      
           
     # Load project tutorial from json file. 
     $ python manage.py api load --json ../biostar-recipes/projects/tutorial.hjson 
     
     # Load recipe jsons in --dir to project --pid. Upload to remote host with -u flag.
     $ python manage.py api load --pid tutorial --dir ../biostar-recipes/recipes/ --recipes
    

## Methods

### Project list 
    GET /api/project/list/

List of projects in a tab delimited fashion with three columns: **ID** , **Name**, and **Privacy**

#### Example
[/api/project/list/](https://www.bioinformatics.recipes/api/project/list/)

    tutorial	Recipe Tutorials	Public
    cookbook	Bioinformatics Cookbook	Public

### Project Information

    GET /api/project/{id}/
    PUT /api/project/{id}/
   
#### Parameters
* _id_: unique project ID

#### Example
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
  

#### Parameters
* _id_: unique project ID

#### Example

[/api/project/image/tutorials/](https://www.bioinformatics.recipes/api/project/image/tutorials/)
    
Image in response:

![20% center](tutorial.png)


### Recipe list
    GET /api/recipe/list/{id}/
    
List of recipes for project in a tab delimited fashion with two columns: **ID** and **Name**

#### Parameters
* _id_: unique project ID

#### Example
[/api/recipe/list/tutorials/](https://www.bioinformatics.recipes/api/recipe/list/tutorials/)

    interface	Interface Elements
    makefile	Makefile Example
    environment	Environment Check

### Recipe Json

    GET /api/recipe/json/{id}/
    PUT /api/recipe/json/{id}/

Recipe JSON used to generate interface

#### Parameters
* _id_: Unique recipe ID

#### Fields in response 
Fields associated with the recipe JSON

#### Example
[/api/recipe/json/starter/](https://www.bioinformatics.recipes/api/recipe/json/starter/)
    
    {
      readlen:
      {
        label: Read Length
        display: INTEGER
        value: 250
        range:
        [
          70
          100000
        ]
      }
      instrument:
      {
        label: Select Instrument
        display: DROPDOWN
        choices:
        [
          [
            hiseq
            Illumina Hiseq
          ]
          [
            pacbio
            Pacific BioSciences Sequel
          ]
          [
            minion
            Oxford Nanopor Minion
          ]
        ]
        value: pacbio
      }
      reference:
      {
        label: Reference Genome
        display: DROPDOWN
        type: FASTA
        source: PROJECT
        value: Genome.fa
      }
      settings:
      {
        name: Starter Recipe
        image: starter.png
        help:
          '''
          This recipe can be a starting point for other recipes.
          # Help
    
          Use this recipe to create new recipes.
          '''
        id: 12
        uid: starter
        project_uid: tutorial
        url: http://localhost:8000
      }
    }

### Recipe Template

    GET /api/recipe/template/{id}/
    PUT /api/recipe/template/{id}/

Recipe template executed during analysis.

#### Parameters
* _id_: Unique recipe ID

#### Example
[/api/recipe/template/starter/](https://www.bioinformatics.recipes/api/recipe/template/starter/)
    
    #
    # A starter recipe with examples.
    #
    
    #
    # You can fill in shell variables
    #
    READLEN=250
    
    echo "Read length: $READLEN"
    
    #
    # Substitute into content
    #
    echo "Referene genome: Genome.fa"
    
    #
    # But you may also use Django Template constructs.
    #
    
    
    echo "Yes, it is Pacific Biosciences!"
    
    #
    # Generate a table of content with all files in the job directory.
    #
    find . -name '*' > files.txt
    
    #
    # Print the contents to the screen
    #
    echo "****** File List: files.txt ****"
    cat files.txt
    # Make a nested directory
    mkdir -p foo/bar
    find . -name '*' > foo/bar/all.txt


### Recipe Image

    GET /api/recipe/image/{id}/
    PUT /api/recipe/image/{id}/
    
    
#### Parameters
* _id_: Unique recipe ID

#### Example

[api/recipe/image/starter/](https://www.bioinformatics.recipes/api/recipe/template/starter/)

Image in response:

![20%](starter.jpeg)