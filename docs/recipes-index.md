# Documentation for Bioinformatics Recipes

## What is the default admin login?

When the site initializes the admin username and password are using the `ADMINS` and the `ADMIN_PASSWORD` settings in `biostar/forum/settings.py`.

By default both the admin login name and the default admin password are set to

    admin@localhost

**Note**: These settings must be changed on a publicly accessible site!

## How to access the Django Admin interface?

* http://127.0.0.1:8000/accounts/admin/

## How to customize the settings?

DO NOT add your custom settings into the public codebase!

The proper practice is to create a separate, independent settings file, then, within that file import **all** default settings. Finally override the fields that you wish to customize in your settings file. For example
create the `my_settings.py` then add into it:

    # Import all default settings.
    from biostar.recipes.settings import *

    # Now override the settings you wish to customize.
    ADMIN_PASSWORD = "foopass"

Apply this settings file with

    python manage.py runserver --settings my_settings.py

Consult the The [Django documentation][django] for details.

[django]: https://www.djangoproject.com/

## How do I deploy the site?

The software follows the recommended practices for developing and deploying [Django web applications][django] .

The [Django documentation][django] contains a wealth of information on the alternative ways to deploy the site on different infrastructure.

Within this setup we recommend the [uwsgi][uwsgi] based deployment.

[uwsgi]:https://uwsgi-docs.readthedocs.io/en/latest/

## How does the site work?

The site is project based. Each project is a collection of data, recipes and results.

Thus each project has three distinct sections:

1. The data.
2. The recipes.
3. The results.

The **Results** are created by applying a **Recipe** on **Data**.

## What is a recipe?

Each recipe is built from two ingredients:

1. The interface specification file.
2. The template specification file.

The **interface** will specify the value of the parameters that get substituted into the **template**.

The **template** contains the commands that need to be executed. The **template** will have
placeholders for the parameter values that the user will need to enter in the interface.

The interface + template will generate a script that the site can execute.

The software will generate an web interface for each parameter specified in the interface. It is this interface where users are able to select the values that their recipe needs to operate.

## Where can I see tutorial recipes?

Visit:

* https://www.bioinformatics.recipes/recipe/list/tutorials/

## Recipe example: Empty Recipe

The simplest recipe is empty for both the **template** and the **data**.

* https://www.bioinformatics.recipes/recipe/view/empty-recipe/

Even though it performs no action it is a valid and working recipe! It demonstrates what takes place when a recipe is run. The results of running the empty recipe are here:

* https://www.bioinformatics.recipes/job/view/a53f6057/

**Note**: To run a recipe you need to have the EXECUTE permission on the project. Admin users automatically have this permission on every project.  If you don't have this permission you
can still see the results that this recipe produces but you would not be able to run the recipe.

Note how even an empty recipe produces a number of outputs. These are files named as follows:

- `run.sh` is the script that executed after being genereated from the template.
- `run.sh.json` contains the data that was used in the template.
- `run.sh.stdout.log` contains the output messages that the recipe produced.
- `run.sh.stderr.log` contains the error messages that the recipe produced.

The log information is also visible on the result page.

### Recipe example: Hello World

Let's develop our recipe a little more.

* https://www.bioinformatics.recipes/recipe/view/hello-world/

In this recipe the template contains the following:

    # This is a regular bash script.

    echo 'Hello World!'

The recipe takes no input and when run simply prints "Hello World". The results of running this recipe can be seen here:

* https://www.bioinformatics.recipes/job/view/3e365b2c/

Note that the words "Hello World" appear on the "Output Messages" tab and are also contained in the file called `stdout.txt`

* https://www.bioinformatics.recipes/job/serve/3e365b2c/runlog/stdout.txt


Make a new recipe and add the following into it:

## The interface

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

#### The template

[template]: https://docs.djangoproject.com/en/1.11/ref/templates/language/

The templates follow the [Django Templating Language Syntax][django] but you would only need
to study that syntax if you need advanced functionality. The most commonly used
functionality that of parameter substitution is  quite straitforward as it substitutes
the value of a parameter if it is enclosed with the `{{ }}` symbols.

Add the following to the template:

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

