# Recipe Tutorial

## Basic Concepts

Each project has three distinct sections:

1. The Data Entries.
2. The Recipes.
3. The Results.

The Results are built by applying a Recipe onto Data Entries.

An interesting functionality of the site is that one Data Entry may be a single 
file or a collection (archive or directory) of any number of files (hundreds, thousands). 
Recipes then can be written to run in parallel on all files in one collection. Hence
a single recipe, in a single invocation can analyze all files in a collection.


For now you don't need to be overly concerned with this feature. 
We just wanted to mention it since it is an interesting and empowering feature.

## What is a Recipe?

Each recipe is built from two ingredients:

1. The interface specification file
2. The template specification file

The **interface**
will specify the value of the parameters that get substituted into the **template**.

The software (`biostar-engine`) generates an interface for each parameter
where the users are able to select the values that their recipe needs to operate.

The **template** contains the commands that need to be executed. The **template** has
placeholders for parameter values that the user will need to enter. 

Go ahead and run the **recipe**. It will be queued and executed shortly.
Investigate the files that it has created.

Let's move onto the second example.

## Recipe: Hello World 1

The simplest recipe is empty for both the **template** and the **data**.

Create it!

The recipe is valid! It just does nothing other then demonstrate
what takes place when an analysis is run. This recipe will be listed as
**No name** and **No summary** since default values will be inserted into
all categories.

Now go ahead and run this analysis.

Note how even an empty recipe produces a number of outputs. There are files named as follows:

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
