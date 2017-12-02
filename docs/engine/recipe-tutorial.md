# Basics

Each project has three distinct sections:

1. The Data Files.
2. The Analysis Recipes.
3. The Analysis Results.

The Analysis Results are built by applying one Analysis Recipe onto a Data File.

One important detail is that a Data File may be one single file or a collection of any number of files (hundreds, thousands).
Many recipes can be run in parallel on all files in one invocation.
For now you don't need to be overly concerned with this feature. We just wanted to mention it since it is such a cool
feature.

# What is a Recipe?

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
Investiage the files that it has created.

Let's move onto the second example.

# Recipe: Hello World 1

The simplest recipe is one that contains **NOTHING** for both the **template** and the **data**.

Even having no content the analysis is valid! It just does nothing other then demonstrate
what takes place when an analysis is run. This recipe will be listed as
**No name** and **No summary** in the [Analysis List](http://www.lvh.me:8000/analysis/list/1/).

Now go ahead and investigate the [**No name**](/analysis/view/4) analysis then
the results of the **No name** analysis.

Follow the [**Show Recipe**](/analysis/recipe/4) action. It shows you how both the **template** and the **data** are empty.
Follow the [**Show results**](/job/list/1/?filter=4) action and click on a result.
Note how even an empty recipe produces outputs. The file named:

- `run.sh` contains the command that were run.
- `run.sh.json` contains the command that were run.
- `run.sh.stdout.log` contains the output messages that the recipe produced
- `run.sh.stderr.log` contains the error messages that the recipe produced

This same information is visible for the job on the previous page.

When a recipe is run a new folder is created that contains a substantial amount of information
on how the recipe worked and what the `biostar-engine` (the software that runs the recipe) deals
with processes.

# Recipe: Hello World 2

Let's make the recipe print **Hello world!**

Visit this recipe and look at the recipe source. Note the changes added to the **template** and the
**data**.

The special entry called **settings** allows you to control
the way the recipe is displayed on the site. This is where the name,
the summary and the description for the recipe comes from.

# Recipe: Hello World 3

This recipe demonstrates how parameters in the **data** get rendered as input widgets.

The tools takes two parameters, organism and cutoff and passes these into the template.
