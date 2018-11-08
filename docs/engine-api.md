# API

### Methods
#### Api list
    GET /api/list

List of recipes with api links corresponding to the JSON and template.

    {
      "R Script":
      {
        uid: 93412cee
        json: /api/recipe/93412cee/json/
        template: /api/recipe/93412cee/template/
      }
      "Makefile Example":
      {
        uid: 16c4f58f
        json: /api/recipe/16c4f58f/json/
        template: /api/recipe/16c4f58f/template/
      }
      "Interface Elements":
      {
        uid: "74391e69"
        json: /api/recipe/74391e69/json/
        template: /api/recipe/74391e69/template/
      }
      "Starter Recipe":
      {
        uid: 4d1846f1
        json: /api/recipe/4d1846f1/json/
        template: /api/recipe/4d1846f1/template/
      }
    }
Only public recipes are shown by default, however provided and api key all recipes are shown.

The Api key is passed as a GET parameter `?k=` in the url, an example: https://www.bioinformatics.recipes/api/list/?k=api-key

### Json

Pick a recipe and click `Json` to view raw JSON of the recipe.

Api links for the JSON look like this: https://www.bioinformatics.recipes/api/recipe/93412cee/json/

Once clicked, the payload will be formatted like such:

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

Pick a recipe and click `Template` to view raw template of the recipe. 

Api links for the template look like this: https://www.bioinformatics.recipes/api/recipe/93412cee/template/

Once clicked, the payload will be formatted like such:

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
