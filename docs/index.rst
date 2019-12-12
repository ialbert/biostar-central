Bioinformatics Recipes: Reproducible Data Analysis
=====================================================================

Bioinformatics Recipes is a `Python <http://www.python.org/>`_ and
`Django <http://www.djangoproject.com/>`_ based data analysis software licensed under the *MIT Open Source License*.
It is a simple, generic, flexible and extensible analytic framework.

Support
-------

The software is open source and free to use under the most permissible license.

The developers of the software are also available to provide commercial level support
for deploying Bioinformatics Recipes for entire organizations. Contact: mailer@bioinformatics.recipes

Requirements: *Python 3.7*

Quick Start
------------

From the biostar source directory::

    # Install the requirements.
    pip install --upgrade -r conf/requirements.txt
    conda install --file conf/conda_requirements.txt

    # Load the environment variables.
    conda activate engine

    # Initialize database, import test data, and run the site.
    make recipes init recipe_demo serve

Visit **http://localhost:8000** to see the site loaded with default settings. Enjoy.

For more information see the documentation below:



.. toctree::
   :caption: Getting Started
   :maxdepth: 2

   recipes/general
   recipes/install

.. toctree::
   :caption: Features
   :maxdepth: 3

   recipes/projects
   recipes/recipes
   recipes/api
   deploy
   customize
   about
   recipes/faq
