# API

## Methods
### Api list
    GET /api/list

List of recipes with api links corresponding to the JSON and template.

#### Fields in response 
JSON dictionary with each recipe keyed by it's `id`.

**key**: Unique recipe ID
* _name_ : Recipe Name
* _json_: API link for the recipe JSON
* _template_: API link for the recipe template

#### Example
[/api/list/](https://www.bioinformatics.recipes/api/list)

    {
      93412cee:
      {
        name: R Script
        json: /api/recipe/93412cee/json/
        template: /api/recipe/93412cee/template/
      }
      16c4f58f:
      {
        name: Makefile Example
        json: /api/recipe/16c4f58f/json/
        template: /api/recipe/16c4f58f/template/
      }
      74391e69:
      {
        name: Interface Elements
        json: /api/recipe/74391e69/json/
        template: /api/recipe/74391e69/template/
      }
      4d1846f1:
      {
        name: Starter Recipe
        json: /api/recipe/4d1846f1/json/
        template: /api/recipe/4d1846f1/template/
      }
    }

### Json

    GET /api/recipe/{id}/json

Recipe JSON used to fill the template

#### Parameters
* _id_: Unique recipe ID

#### Fields in response 
Fields associated with the recipe JSON

#### Example
[/api/recipe/93412cee/json](https://www.bioinformatics.recipes/api/recipe/93412cee/json/)

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

    GET /api/recipe/{id}/template

Recipe template executed.

#### Parameters
* _id_: Unique recipe ID

#### Fields in response 
Plain text response of the recipe template.


#### Example
[/api/recipe/93412cee/template](https://www.bioinformatics.recipes/api/recipe/93412cee/template/)

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
