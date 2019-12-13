# FAQ - Frequently Asked Questions 

## Answers to reviewers

This section was created specifically to answer questions reviewers posed while reviewing our scientific publication.

### Who can run a recipe?

To run a recipe the recipe must be **authorized**  and the user must have **trusted** designation.

### What is an authorized recipe?

Since a recipe may contain shell commands or other code security checks are needed to avoid the misuse of computational resources. For that reason every new recipe starts out in a so-called **pending authorization** state, labeled with an orange ribbon. These recipes cannot be executed on the website, but may be inspected, viewed, shared or downloaded.

![](images/authorization-pending.png)

A user with administrative privileges (an administrator) must approve a recipe (in the recipe edit window) for the recipe to become executable within the website. A green ribbon decorates authorized recipes.

![](images/authorization-valid.png)

### Who is a trusted user?

Each user has a designation: **trusted** or **visitor** that controls their ability to run recipes. Even when it comes to **approved** recipes, only users with **trusted** designation may run these recipes.

![](images/authorization-untrusted-user.png)

![](images/user-all-menu-items.png)
![](images/user-admin-icon.png)


### How does one become a trusted user?

The site administrators are able to change the designation of any user.

The restrictions that we have in place provide high granularity control of the computational resources.

The owner of the site decides which users and which recipes gain the privileges to use the computational resources. Other groups running the Recipes software may set up their system in such a way to automatically trust every new user that signs up, and they may also choose to approve every recipe that is created automatically.

### Can this software be run on my machine?

The software was designed with decentralization in mind. The software runs on any operating system: Linux, macOS and Windows 10 (with Linux Subsystem), and on any hardware that supports Python. We routinel run the software on a MacBook Air laptop, on a single computer serving a lab, and on a high-performance multicore server.

Installation takes little more than minutes and requires no special software, just the ability to run the Python programming language.  We do envision different groups running their personalized instances of the software to serve local needs.

### Can the software be used in bioinformatics education?

Even though originally the platform was designed to provide bioinformatics support to our biologist collaborators we have also found that the use of recipes integrates exceedingly well within bioinformatics curricula.

Specifically, when covering more advanced topics, educators must present a series of commands that produce data analysis steps. Yet currently, there is no straightforward way to publish both the code and the results in one location. We can use say GitHub to publish code, but we can't use GitHub to execute the code, and we can't use GitHub to store large datasets either.

 In contrast, when using bioinformatics recipes, a student can see both the code and all the files and results generated while running the code. In addition, they can also view results generated with different runtime parameters; all results are linked to the recipe that produced them.

Finally, students can readily copy the recipe over to their projects, make changes to it and see the results of their changes, all in the same interface.

### The software seems conceptually most similar to Galaxy, but how does it differ?

When compared to Galaxy the Recipes framework presents several fundamental differences in its operating principles:

1. The framework was designed to capture entire data analysis pipelines. The Galaxy framework was designed to represent individual software tools.
1. The code downloaded from the the site can be run on any other computer (that has the software installed) and will produce the same output as seen via the web interface.
1. A Galaxy analysis can only be run inside the Galaxy software. Recipes are designed to be shared and expanded upon by various users.
1. Users may create different interfaces for recipes copied from someone else. The end-user cannot change interfaces in Galaxy.
1. The framework is also a data analysis know-how, social interaction and training material distribution framework.

### What are the advantages of the recipes over Galaxy?

One major advantage, in our opinion, is the independence of the method from the platform.

A bioinformatics recipe is an independent piece of code that educates and trains bioinformatician how to perform the analysis themselves. The recipe may be run either on the recipes website or on any other  computational infrastructure. The code obtained are recipes are identical to code a bioinformatician would develop at the command line.

Additionally the platform also serves as knowledge distribution. Users have the ability to build upon each others' know-how and expertise. A user may take an existing multistage analysis and add/remove/customize that analysis for their needs.

### If you have a local Galaxy server and a local Recipe instance - what functional differences do users see?

It is possible to set up recipes that look and behave like tools in Galaxy. For example, if one were to wrap individual tools into recipes, then from the usability perspective, Galaxy and Recipes would be nearly identical.

In that sense, the Recipes website is a superset of some of Galaxy's functionality.
That being said, this is not how we envision using the Recipes, we advocate building pipelines rather than individual tools, a model that is not directly applicable in Galaxy.

### How does the performance compare?

Our frameworks uses Python 3 and was developed with [Django][django], a well documented platform with extensive use in the information technology industry.

When it comes to the performance of the analyses themselves, these depend on the choice of methods and on the infrastructure that runs the server.

[django]: https://www.djangoproject.com/

### This application is designed to serve individual groups, but groups on what scale?

Django, the application server that our application is built upon has wide acceptance in the information technology industry. Django runs platforms such as Pinterest and Instagram.  Thus we believe that software can be made to scale up to support computing at a supercomputer scale.

Our current focus, based on the priorities of the funding that we have received, was to develop a system that serves groups consisting of dozens of bioinformaticians interacting with hundreds of end-users.

### Can it cope with multiple users, for example, what are the minimum requirements for installing the web application?

The bioinformatics recipes software itself has extremely-low memory and CPU overhead.
We estimate a few hundred megabytes and less than 1% CPU utilization.

Of course, when we run an analysis, the resource utilization depends on the tasks in the processes that are employed - what is important to note is that our software imposes minimal overhead. There are no constraints for the number of users that can access/read/share/ copy recipes. Hundreds of thousands of users could be browsing the recipes with millions of page views every month.

### Where (and how) recipes actually run when you execute them through the web application?

The recipes are currently executed on the same platform that runs the webserver. Since the service itself has minimal overhead, the entire computational infrastructure is available for computation.

The recipes software comes with a built-in job queuing system that will honor this setting.
Thus site owners can set up more or fewer simultaneous job queues depending on their computational resources.

### How are new recipes added to the system?

Recipes are created in two ways: either as new recipes by selecting the "Create recipe" button or by "copying" or "cloning" an existing recipe (we provide the explanation for the two terms in the next answer).


### If different people are developing recipes, how do you standardise the way the recipes are created?

This problem is both important yet somewhat of a challenge to implement in a way that is not overbearing and limiting yet maintains utility to the users. 
We have chosen the model described below, but we are open to suggestions from our community and may revisit the implementation later. 
For that reason, we would rather maintain this information in the main documentation rather than in the publication itself as it is possible that later on, our implementation changes.

### Every recipe that a user may access can be duplicated in a new project as a "clone" or as a "copy".

A cloned recipe remains in sync with the original recipe that it was cloned from. 
Clones cannot be changed and track the original recipe. A change to the original recipe will immediately be reflected in all the clones. The purpose of a cloned recipe is to ensure that a recipe is the same across multiple projects and individuals.

![](images/paste-as-clone.png)

The second method to duplicate a recipe is to copy it. A copied recipe is a brand new recipe filled with the content from an existing recipe. 
When a recipe is copied the provenance to the original recipe is not maintained. It becomes the responsibility of the author of the recipe to maintain the relevant information in the documentation of the recipe. 

![](images/paste-as-new.png)


### What conventions should be followed, how should they be documented, what is the minimum amount of provenance that the scripts should produce, and how should that be presented to the users?

Choosing the appropriate levels of documentation and provenance is a difficult question that we still debate and discuss. 
We would like to avoid being either too lax or too stringent. The concepts that we popularize in this tool are new; the approach is different from past models. 
In the paragraphs above, we discuss many of the challenges that we face. 

We do plan to evolve our views as needed. Currently, we chose to approve only the recipes where the documentation is appropriate, and provenance is properly noted. 
Hence the "approved" state of a recipe is a manually curated process. We hope that with time, as we better understand the standards and requirements, we will be able to automate the process.

## Other questions
