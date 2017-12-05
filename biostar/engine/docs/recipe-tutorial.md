# Basic Concepts

Each project has three distinct sections:

1. The Data.
2. The Recipes.
3. The Results.

The **Results** are created by applying a **Recipe** on **Data**.

## What is a Recipe?

Each recipe is built from two ingredients:

1. The interface specification file.
2. The template specification file.

The **interface** will specify the value of the parameters that get substituted into the **template**.

The **template** contains the commands that need to be executed. The **template** will have
placeholders for the parameter values that the user will need to enter in the interface.

The interface + template will generate a script that the site can execute.

The **Biostar Engine** will generate an web interface for each parameter specified in the interface.
It is this interface where users are able to select the values that their recipe needs to operate.

## Recipe: Hello World 1

The simplest recipe is empty for both the **template** and the **data**.

Even though it is empty it is a valid and working recipe! It just does nothing other then demonstrate
what takes place when a recipe is run. 

To run a recipe you need to have the EXECUTE permission. If you don't have this permission you 
can still see the results that this recipe produces by clicking the `Go to Results` link
that will filter the results to show you only those that were created by this recipe.

Let's investigate one result.

Note how even an empty recipe produces a number of outputs. These are files named as follows:

- `run.sh` is the script that executed after being genereated from the template. 
- `run.sh.json` contains the data that was used in the template.
- `run.sh.stdout.log` contains the output messages that the recipe produced.
- `run.sh.stderr.log` contains the error messages that the recipe produced.

The log information is also visible on the result page.

## Recipe: Hello World 2

Let's add something more to each field.

Make a new recipe and add the following into it:

### Interface

Let's add some parameters to the recipe. Enter the following:

    {
        foo: {
            value: 100
        }
      
    bar: {
            value: 200
        }
    }

[hjson]: https://hjson.org/
[json]: https://en.wikipedia.org/wiki/JSON

The syntax follows the so called JSON notation. We are using
a variant of JSON that is better suited for human input
called [HJSON][hjson] (Human JSON). The HJSON variant 
is an extension of [JSON][json] that is fully compatible 
with JSON so you may use the original [JSON][json] notation 
if you so desire.
 

All data objects are dictionaries. The outer (root) dictionary `{}` is keyed with `foo` and `bar`.
In our nomenclature `foo` and `bar` will be the *parameters* to the template. 
Each parameter holds another dictionary `{}` and parameter has  a `value` key in this dictionary.

### Template

[template]: https://docs.djangoproject.com/en/1.11/ref/templates/language/

The templates follow the [Django Templating Language Syntax][django] but you don't really need
to study that syntax unless you need some advanced functionality. The most commonly used
functionality is self explanatory.

Now add the following to the template:

    echo FOO={{foo.value}}
    echo BAR={{bar.value}}

Press the **Preview** button and your generated script then will produce

     echo FOO=100
     echo BAR=200
    
Note how our interface does not have entires for `foo` and `bar`. This 
is because it does not understand how to generate a widget for `foo` and `bar`. In the 
next tutorial we will show you that.

For now go ahead, **Save** the recipe then go and run it. Once the run 
completes evaluate the results. What do you see?

## Recipe: Hello World 3

This recipe demonstrates how parameters in the **data** get rendered as input widgets.

The tools takes two parameters, organism and cutoff and passes these into the template.

## Recipe: Hello World 4

Let's make the recipe print **Hello world!**

Visit this recipe and look at the recipe source. Note the changes added to the **template** and the
**data**.

The special entry called **settings** allows you to control
the way the recipe is displayed on the site. This is where the name,
the summary and the description for the recipe comes from.
