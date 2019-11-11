## Documentation for Bioinformatics Recipes

### What is the default admin login?

When the site initializes the admin username and password are using the `ADMINS` and the `ADMIN_PASSWORD` settings in `biostar/acccounts/settings.py`.
 
 By default both the admin login name and the default admin password are set to

    admin@localhost

**Note**: These settings must be changed on a publicly accessible site!

### How to access the Django Admin interface?

* http://127.0.0.1:8000/accounts/admin/

### How to customize the settings?

DO NOT add your custom settings into the public codebase!

The proper practice is to create a separate, independent settings file, then, within that file import **all** default settings. Finally override the fields that you wish to customize in your settings file. For example
create the `my_settings.py` then add into it:

    # Import all default settings.
    from biostar.recipes.settings import *

    # Now override the settings you wish to customize.
    ADMIN_PASSWORD = "foopass"

Apply this settings file with

    python manage.py runserver --settings my_settings.py

Consult the [Django documentation][django] for details.

[django]: https://www.djangoproject.com/

### How do I deploy the site?

The software follows the recommended practices for developing and deploying [Django web applications][django] .

The [Django documentation][django] contains a wealth of information on the alternative ways to deploy the site on different infrastructure.

Within this setup we recommend the [uwsgi][uwsgi] based deployment.

[uwsgi]:https://uwsgi-docs.readthedocs.io/en/latest/

### How does the site work?

The site is project based. Each project is a collection of data, recipes and results.

Thus each project has three distinct sections:

1. The data.
2. The recipes.
3. The results.

The **Results** are created by applying a **Recipe** on **Data**.

### What is a recipe?

Each recipe is built from two ingredients:

1. The interface specification file.
2. The template specification file.

The **interface** will specify the value of the parameters that get substituted into the **template**.

The **template** contains the commands that need to be executed. The **template** will have
placeholders for the parameter values that the user will need to enter in the interface.

The interface + template will generate a script that the site can execute.

The software will generate an web interface for each parameter specified in the interface. It is this interface where users are able to select the values that their recipe needs to operate.

### Where can I see tutorial recipes?

See the url below for a number of recipes of increasing complexity:

* https://www.bioinformatics.recipes/recipe/list/tutorials/

### Recipe example: Empty Recipe

The simplest recipe is empty for both the **template** and the **data**.

* https://www.bioinformatics.recipes/recipe/view/empty-recipe/

Even though it performs no action it is a valid and working recipe. Its purpose is to demonstrate what takes place when a recipe is run. The results of running the empty recipe are here:

* https://www.bioinformatics.recipes/job/view/a53f6057/

**Note: You need to be a trusted user to run a recipe**. Admin users automatically have this permission on every project.  If you don't have this permission you
can still see the results that this recipe produces but you would not be able to run the recipe.

Note how even an empty recipe produces outputs. These are files named as follows:

- `recipe.sh` file is the script that executed after being generated from the template.
- `runlog/input.json` file contains the data that was used in the template.
- `runlog/stdout.txt` file contains the output messages that the recipe produced.
- `runlog/stderr.txt` file contains the error messages that the recipe produced.

The contents of `stdout.txt` and `stderr.txt` are also visible on the result page.

### Recipe example: Hello World

Let's write a recipe that prints "Hello World" to the screen.

* https://www.bioinformatics.recipes/recipe/view/hello-world/

In this recipe the template contains the following:

    # This is a regular bash script.

    echo 'Hello World!'

The recipe is a bash script that prints "Hello World" to the screen. The results of running this recipe can be seen here:

* https://www.bioinformatics.recipes/job/view/3e365b2c/

Note that the words "Hello World" also appear on the "Output Messages" tab and are contained in the file called `stdout.txt`

* https://www.bioinformatics.recipes/job/serve/3e365b2c/runlog/stdout.txt

Make a new recipe and add the following into it:

### Recipe example: Download FASTQ data by SRA number

Suppose we wish to create a recipe that downloads and unpacks FASTQ data from the short read archive.
The code we wish to deploy is:

    # The SRR run number.
    SRA=SRR519926

    # Download 1000 reads from SRA.
    fastq-dump --split-files -X 1000 $SRA

but we want to make the selection of the SRA number controllable by the user.

We start by copying over any other existing recipe. Start with the "empty recipe" for example.

Find the "Interface link"  it is in `More -> Interface` then paste the code above into the template section. Click "Preview" to see what the code will look like, in this case since the code does not have any modifiable region it will look the same after the preview.

Save this recipe. You have recipe that works on one specific SRA number. If that is all you wanted you would be done with the recipe.

To make the input overrideable we need to add the following to the Interface JSON section (this might be already filled out to some default settings. Replace all that with:

    {
        settings: {

        }

        sra: {
            value: SRR519926
        }
    }

All data objects are dictionaries. The `settings` key is internal. The `sra` key is a parameter to the script. To access this parameter from the script change the template to

    # The SRR run number.
    SRA={{ sra.value }}

    # Download 1000 reads from SRA.
    fastq-dump --split-files -X 1000 $SRA

Note here that we access the value of the parameter `sra` with ``{{sra.value}}``.

If you preview your recipe again you will see that it produces the same output as before. The value is filled into the script automatically.

But the interface is still empty as the site does not yet know how to render a graphical widget to the parameter. To tell the site how to render the parameter expand the interface JSON to look like this:


    {
        settings: {

        }

        sra: {
            display: TEXTBOX
            value: SRR519926
            help: An SRA run number
            regex: \w{1,9}$
        }
    }

When you press the "Preview" again you will see the following interface:


<img src="https://raw.githubusercontent.com/ialbert/biostar-central/master/docs/recipes/interface-1.png" width="600">

With have instructed the site to display the parameter as a `TEXTBOX` that only accepts a single maximum 9 letter word as input. It also renders a small help under the textbox to inform the user of the purpose of the input.

And that's it! Save the recipe, and now you have just written a simple recipe that others may run and reuse.

<img src="https://raw.githubusercontent.com/ialbert/biostar-central/master/docs/recipes/interface-2.png" width="600">


### What format is the interface in?

The JSON syntax follows  a variant of JSON that is better suited for human input
called [HJSON][hjson] (Human JSON). HJSON
is an extension of [JSON][json] that is fully compatible
with JSON so you may use the original [JSON][json] notation
if you so desire.

### Where can I see more code examples for interface and scripts?

Visit the recipes website and see the various example recipes:

* https://www.bioinformatics.recipes/
