
Introduction
-------------

The codebase requires Python_ (version 2.6 or higher) and Django_ (version 1.3 or higher) to run.

We also use JQuery_ for javascript and `markitup`_ for 
rich text javascript editor. These do not need to be installed and will be loaded over the internet.

Installation
------------

To install Django 1.3 please see the INSTALL file.

There is a main run manager in the root directory::

    $ biostar.sh 

execute it with no parameters for information on usage. This run manager 
can take one or ore commands. For example to initialize, populate and run the server
one would invoke it as follows::

    $ biostar.sh init 
    $ biostar.sh populate
    $ biostar.sh run

Alternatively one may run all these commands all at once::

    $ biostar.sh init populate run

**Note**: If the models change you must reset and reinitialize your database::

    $biostar.sh delete init populate

The default server will bind the all IP adapters (0.0.0.0) and port 8080. Visit http://localhost:8080 to see
interact with your version of the test server. Edit the `biostar.sh` script to override the various settings.

Layout
------

The Python code is found in the biostar folder. Templates, static content 
(css, images, javascript) and data are found in the home directory. 
There is partial datadump of the existing BioStar content in the 
`home/import/datadump`. The `populate` command will load 
this data into the current database.

Colorscheme
-----------

  * Purple: `#8F2C47`
  * Green: `#75845C`

.. _Django: http://www.djangoproject.com/
.. _Python: http://www.python.org/
.. _JQuery: http://jquery.com/
.. _markitup: http://markitup.jaysalvat.com/home/
