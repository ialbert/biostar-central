# Social authentication

&bull; back to [Documentation Home](index.md)

The social logins settings will need to be initialized with the proper authentication
parameters. In most cases this involves setting up an account at the provider then
obtaining a ``CLIENT ID`` and a ``CLIENT SECRET``.

Once the credential information is known Biostar users with staff level permissions can
access the Django admin at ``http://yoursite.here/admin`` and add these credentials
into the ``Social Applications`` table.

Note: The authentication information is stored in the database.
Database dumps will contain this potentially sensitive information.

To avoid having to manually add social account settings users can export social
authentication related values from their database with:

    python manage.py dumpdata socialaccount.socialapp > socialapp.json

Loading these settings would be performed with:

	python manage.py loaddata socialapp.json

# Setup Examples

Biostar uses ``django-allauth`` for authentication. The documentation for it can be found at
http://django-allauth.readthedocs.org/en/latest/providers.html

We demonstrate a few simple examples below

## Adding Google authentication:

* Ensure that ``allauth.socialaccount.providers.google`` is listed in your ``INSTALLED_APPS``
* Access the Google Developer Console: https://cloud.google.com/console/project
* Find the section called ``API&auth``
* Create a new Client ID.
* Ensure that the ``REDIRECT URIS`` is set to ``http://yoursite.here/accounts/google/login/callback/``
* Add the ``CLIENT ID`` and ``CLIENT SECRET`` into your ``Google Social Account`` on the Admin site.

## Adding Persona Authentication:

Setting up persona authentication is somewhat different from others.

* Ensure that ``allauth.socialaccount.providers.persona`` is listed in your ``INSTALLED_APPS``
* Ensure that an audience is set up in your settings:

		SOCIALACCOUNT_PROVIDERS = {
			'persona': {
				'AUDIENCE': 'http://yoursite.here/',
				'REQUEST_PARAMETERS': {'siteName': 'Your Site Name'}
			}
		}
