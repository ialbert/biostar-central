# Deploying Biostar

&bull; back to [Documentation Home](index.md)

All settings are specified in the ``biostar3/settings/values.py`` and ``biostar3/settings/base.py`` files.

Ensure that your settings file inherit from the ``biostar3.settings.base`` file. This will also load the ``biostar3/settings/values.py``
file. A typical settings file will start like so

    # Import all values from the base then override site specific settings.
    from biostar3.settings.base import *

Then server admins are expected to selectively override settings as needed.

There are example settings in the the ``run/*`` folder for example see the ``run/sqlite.py`` or ``run/postgres.py`` 
file.

The important settings that must be overridden for the site to work are:
	
	# Set up the default site
	SITE_ID = 1
	SITE_NAME = "Biostars Q&A"
	SITE_DOMAIN = "www.lvh.me:8080"
	SITE_SCHEME = "http"
	
	# Set up the moderator site.
	MODERATOR_SITE_NAME = "Biostars Moderators"
	MODERATORS_SITE_DOMAIN = "moderators.lvh.me:8080"
	
	
