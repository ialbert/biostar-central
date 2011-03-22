
Introduction
-------------

The codebase requires python 2.6 or higher and Django version 1.2 or higher.
Using the Blueprint CSS framework for styling. 

There is a main run manager in the root directory:

    $ biostar.sh 

execute it with no parameters for information on usage. A typical initial run would be:

    $ biostar.sh init populate run

Layout
-------

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