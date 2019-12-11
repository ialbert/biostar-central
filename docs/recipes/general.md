# General Concepts

## What is a recipe?

Each recipe is built from two ingredients:

1. The interface specification file.
2. The template specification file.

The **interface** will specify the value of the parameters that get substituted into the **template**.

The **template** contains the commands that need to be executed. The **template** will have
placeholders for the parameter values that the user will need to enter in the interface.

The interface + template will generate a script that the site can execute.

The software will generate an web interface for each parameter specified in the interface. It is this interface where users are able to select the values that their recipe needs to operate.

A recipe consists of a "JSON definition file" and a "script template".

The simplest JSON definition file is

    {}

A simple script template can be completely empty or might contain just:

    echo 'Hello World!'
    

## Who can execute recipes?


Only **admins**,**staff**,and **trusted users** can execute recipes. 

Trusted users are picked by admins.


## How are recipes standardized?

(TODO)
 
All recipes are built from two ingredients: 

1. The interface specification file.
2. The template specification file.

Standardizing involves the putting of the samples in the  


If different people are developing recipes, 

how do you standardise the way the recipes are created? 

## Where do recipes run when executed?

Where do recipes actually run 


## Conventions when creating recipes


## Documenting recipes


## Amount of provenance that recipes produce


## Presenting the recipe provenance and result.


 How does Bioinformatics Recipes differ from Galaxy


 Advantages 


  Functional differences seen by users


 Performance comparison 


 Scaling to larger users




