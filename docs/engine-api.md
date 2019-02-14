# API

## Methods

### Project list 
    GET /api/project/list/

List of projects in a tab delimited fashion with three columns: **ID** , **Name**, and **Privacy**

##### Example
[/api/project/list/](https://www.bioinformatics.recipes/api/project/list/)

    tutorial	Recipe Tutorials	Public
    cookbook	Bioinformatics Cookbook	Public


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

#### Parameters
* _id_: Unique recipe ID

#### Fields in response 
Fields associated with the recipe JSON

#### Example
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

### Template

    GET /api/recipe/{id}/template/
    PUT /api/recipe/{id}/template/

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
