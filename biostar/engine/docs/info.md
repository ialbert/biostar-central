## What is a recipe?

A recipe is a data analysis script. It may be a `bash` script, an `R script`, a `Makefile`, 
or any other command line tool. 

## Do recipes represent a single software tool?

No. 

A recipe may represent a single software tool and it could also represent an entire *data analysis pipeline*.

## Can recipes be run outside of the site?

Yes. 

A recipe is a command line construct that could be run just 
as well from the command line on any other computer.

Recipes does not favor any particular approach or methodology and can be integrated 
with a wide variety of other tools and technologies. 
   
## How do I use recipes?

You can copy a recipe to a new project that you have created. Once copied the 
recipe it becomes "yours". You can do anything with it without affecting the original recipe.

## How do I create results?
  
When a recipe is run a result row is created. Each result represents a directory that in turn
will contain the of files that the recipe has created while running.

# Are there different access rights?

Yes. 

Users on the site can have different roles. 
Depending on the role some actions may not be permitted.

Project owners may add people to their projects with different roles.

## What kinds of access roles exist?

Project owners have the option to limit access to their project at the following
levels of granularity:

1. No access to project (private project)
2. Read project content
3. Read recipe code in project
4. Execute recipe in project
5. Edit information in project
6. Upload data into project
7. Administer project

The permissions are hierarchical. Each rank includes all permissions at lower rank. 
For example the `Upload data` permission will also include `Execute Recipe` permissions (and all other
lower ranks).

## May I modify a recipe?

You can adapt and remix recipes to create new and modified ones.

## Can I run a modified recipe?

Since the site runs you code unaltered a modified recipe must be reviewed and approved by 
a site administrator before you can run it yourself.

## Can I run a modified recipe without approval?

Modifications by users that have a site administrative role
are approved automatically.

## How can I become a site administrator?

You may run your own version of the Biostar Engine and you may give yourself 
administrative role on your own site. 
   
## Where do I find more information about the software?

The software that runs this site is called the **Biostar Engine**
and can be obtained from the [biostar-engine][engine] repository . 

[engine]: https://github.com/biostars/biostar-engine
[recipes]: https://github.com/biostars/biostar-recipes

## How do I make new recipes?

Visit the [biostar-recipes][recipes]
repository for information on how to create new recipes.

[recipes]: https://github.com/biostars/biostar-recipes
  
# What is a Public project?

The *Public* project designation allows `Read recipe` access to every visitor to the site.

# What is a Private project?

The *Private* project designation limits access to people 
added to the project by the project administrator.
  
# How do I upload large files?

A secure FTP server integrated with the site allows the upload and download of large datasets.

