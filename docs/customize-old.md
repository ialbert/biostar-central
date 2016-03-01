## Customize

**Note:** when customizing any part of Biostar users do not typically need to modify the original files. 

Typically one would create a new folder or settings file and instruct Biostar to load those instead of the origina file. This way changes are easier to pull.

### Customizing settings
[django]: https://www.djangoproject.com/

Biostar is a [Django][django] based application and follows the Django conventions for loading settings. Upon starting up the value of the `DJANGO_SETTINGS_MODULE` environment variable determines which settings module will be imported.

To more easily change these settings the default Biostar configuration will also read more environment variables.


- - -

The section below is changing

- - - 
 
### Getting started


To customize Biostar you will need to overwrite the settings.

The settings are located in two files, an environment file and a python module.

Technically speaking only the python module is absolutely required. We chose to
put certain variables into an environment so that it is easier to share the same
values across multiple settings as well as to have them available at the command line.

The capitalized variables defined in the ``biostar/settings/base.py`` file directly
affect the operation of the site. Each group of variables is fairly well documented
to describe what it is used for. Typically only a small subset will ever need
to be changed.

A typical run sources a shell program and loads variables. The python module then looks
for certain variables in the environment and uses those to set its parameters.

Custom modules
--------------

To get started create a new empty python file. Say ``custom.py``. This file will govern
the entire operation of your site and is the so called  **django settings module**.

Place your django settings module into a folder that python will recognize as a
package directory (has a ``__init__.py`` in it).
For example we use the ``live`` directory in the biostar source. The file will then be
located in ``live/custom.py``. Into this new django settings module place the following::

    from biostar.settings.base import *

This line will ensure that all the default variable are loaded from the base module. We
can now selectively override one or more of these base variables.

By default the emails are printed to the console.
To use a typical SMTP based email service add the following to your custom settings module so
that the entire file looks like::

    from biostar.settings.base import *
    EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'

We are done with django settings module. Now an environment file needs to be created.

Copy the ``conf/defaults.env`` file to a new file. The name typically matches
the settings file that you want to activate. Let's call it ``live/custom.env``::

    cp conf/defaults.env live/custom.env

Modify the ``custom.env`` file so that in it the ``DJANGO_SETTINGS_MODULE``
variable points to your custom django settings module ``live/custom.py``.
This needs to follow the Python convention for import, this means to use a dot instead of
a slash like so ``live.custom`` find and override the line to look like this::

    export DJANGO_SETTINGS_MODULE=live.custom

In this environment file you may also override other variables. Notably
find the ``EMAIL_HOST``, ``EMAIL_PORT``, ``EMAIL_USER``, ``EMAIL_PASSWORD`` variables and
fill in the values that are specific to your internet provider. If you get an
error sending emails then this information is not set properly.

To use the new environment to start
the site by sourcing this script instead of the default one.
You will need to run this once per terminal::

    source live/custom.env

To test that the email was set up correctly::

    python manage.py test_email

This command will send a test email to the email address listed in ``ADMIN_EMAIL`` in the environment file.

A typical site initialization would be::

    source live/custome.env
    ./biostar.sh delete init import run

There are no limitations how many settings one may have. To check which environment is loaded run the
``biostar.sh`` manager on its own. The last printout will display the current django settings module.

Advanced example
----------------

To change the
upper navigation bar ook inside ``biostar.settings.base`` find and copy the following
categories into your custom django settings module ``custom.py`` file then modify as wish::

    # The start categories. These tags have special meaning internally.
    START_CATEGORIES = [
        "Latest", "Tutorials", "Tools",  "Jobs", "Forum", "Unanswered",
    ]

    # These should be the most frequent (or special) tags on the site.
    NAVBAR_TAGS = [
        "Assembly", "RNA-Seq", "ChIP-Seq", "SNP",
    ]

    # The last categories. Right now are empty.
    END_CATEGORIES = [

    ]

    # These are the tags that will show up in the tag recommendation dropdown.
    POST_TAG_LIST = NAVBAR_TAGS + ["software error"]

    # This will form the navbar
    CATEGORIES = START_CATEGORIES + NAVBAR_TAGS + END_CATEGORIES

See the ``conf/defaults.env`` and for all the parameters that can to be customized.

The ``SITE_STYLE_CSS`` and ``SITE_LOGO`` settings permit loading up custom sytles.
See the ``/static/themes`` folder
for examples.



