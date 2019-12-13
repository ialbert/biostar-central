Bioinformatics Recipes: Reproducible Data Analysis
=====================================================================

Bioinformatics Recipes is a `Python <http://www.python.org/>`_ and `Django <http://www.djangoproject.com/>`_ based data analysis software licensed under the *MIT Open Source License*.
It is a simple, generic, flexible and extensible analytic framework.

Quick Start
------------

From the biostar source directory::

    # Install the requirements.
    pip install -r conf/requirements.txt

    # Initialize database, import test data, and run the site.
    make recipes init demo serve

Visit **http://localhost:8000** to see the site loaded with default settings. Enjoy.


For more information see the documentation below:

.. toctree::
   :caption: Getting Started
   :maxdepth: 1
   :hidden:

   recipes/faq
   recipes/install

.. toctree::
   :caption: Features
   :maxdepth: 1
   :hidden:

   recipes/projects
   recipes/recipes
   recipes/api
   deploy
   customize
   about

