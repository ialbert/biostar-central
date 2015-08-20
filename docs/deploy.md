# Deploying Biostar

&bull; back to [Documentation Home](index.md)

All settings are specified in the ``biostar3/settings/values.py`` and ``biostar3/settings/base.py`` files.

Ensure that your settings file inherit from the ``biostar3.settings.base`` file. This will also load the 

    # Import all values from the base then override site specific settings.
    from biostar3.settings.base import *

There are example settings in the the ``run/*`` folder for example see the ``run/sqlite.py`` or ``run/postgres.py`` 
file.

Then selectively override 