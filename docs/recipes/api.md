# API

## Commands 

Pull API:

    $ python manage.py api pull --help
    
    optional arguments:
      -h, --help     show this help message and exit
      -r, --recipes  Pull recipes of --pid
      --url URL      Site url.
      --key KEY      API key. Required to access private projects.
      --rid RID      Recipe uid to dump.
      --pid PID      Project uid to dump.
      --dir DIR      Directory to store in.

    # Dump project from remote url (--url).
    $ python manage.py api pull --pid tutorial --dir tmp/remote/projects/ --url URL

    # Dump recipes from remote url ( --url) 
    $ python manage.py api pull --pid tutorial --dir tmp/remote/recipes/ --url URL --recipes
    
    
Push API:
    
    $ python manage.py api push --help
    optional arguments:
      -h, --help           show this help message and exit
      -u, --url_from_json  Extract url from conf file instead of --url.
      --url URL            Site url.
      --key KEY            API key. Required to access private projects.
      --rid RID            Recipe uid to load.
      --pid PID            Project uid to load.
      --dir DIR            Directory with json files to load in bulk.
      --json JSON          Project or recipe JSON file to load.
                      
           
     # Load project tutorial from json file. 
     $ python manage.py api push --json ../biostar-recipes/projects/tutorial.hjson 
     
     # Load recipe jsons in --dir to project --pid. Upload to remote host with -u flag.
     $ python manage.py api push --pid tutorial --dir ../biostar-recipes/recipes/ -u --key API_KEY
    

## Methods

### Listing 
    GET /api/list/

List projects and recipes in a tab delimited fashion with columns: **Project ID** , **Project Name**, **Recipe ID**, **Recipe Name**, **Privacy**

#### Example
[/api/list/](https://www.bioinformatics.recipes/api/list/)

    tutorial	Recipe Tutorials	environment	Environment Check	Public
    tutorial	Recipe Tutorials	interface	Interface Elements	Public
    tutorial	Recipe Tutorials	makefile	Makefile Example	Public
    tutorial	Recipe Tutorials	rscript	R Script	Public
    cookbook	Bioinformatics Cookbook	quality-check	Improve the quality of sequencing reads	Public
    cookbook	Bioinformatics Cookbook	fastqc	Visualize FASTQ data quality	Public
    cookbook	Bioinformatics Cookbook	pseudo-alignment	RNA-Seq Transcript Abundance Estimation	Public
    cookbook	Bioinformatics Cookbook	augustus	Gene Prediction	Public


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

![20% center](images/tutorial.png)


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

![](images/starter.jpeg)



### Data Update

Adding data api docs.

` curl  --form “@file=/Users/natay/Desktop/apps/biostar-central/SimpleWorkflowMNIST.ipynb”  
        http://localhost:8000/api/data/data-1/`
