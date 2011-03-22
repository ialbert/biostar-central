
Introduction
-------------

The codebase requires Python_ (version 2.6 or higher) and Django_ (version 1.2 or higher).
Currently using the `Blueprint CSS`_ framework for layout and JQuery_ for javascript. 

Installation
------------

First unpack the import data located in `home/import`::

    $ cd home/import
    $ unzip -d output biostar-20100415142324.zip

There is a main run manager in the root directory::

    $ biostar.sh 

execute it with no parameters for information on usage. A typical initial run would be::

    $ biostar.sh init populate run

Layout
------

The Python code is found in the biostar folder. Templates, static content 
(css, images, javascript) and data are found in the home directory. 
There is partial datadump of the existing BioStar content in the 
home/import/ directory. Unzip it before populating the data. 
The populate command will load this data into the current database.

Notes
-----

Upon the first run a secret-key.txt file is created in the home directory. 
The content of this file will serve as default administrative password.

Colorscheme
-----------

Purple: #8F2C47
Green: #75845C

.. _Blueprint CSS: http://www.blueprintcss.org/
.. _Django: http://www.djangoproject.com/
.. _Python: http://www.python.org/
.. _JQuery: http://jquery.com/