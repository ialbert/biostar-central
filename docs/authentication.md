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

# Examples

## Adding Google authentication:

* Ensure that ``allauth.socialaccount.providers.google`` is listed in your ``INSTALLED_APPS``
* Access the Google Developer Console: https://cloud.google.com/console/project
* Find the section called ``API&auth``
* Create a new Client ID.
* Ensure that the ``REDIRECT URIS`` is set to ``http://yoursite.here/accounts/google/login/callback/``
* Add the ``CLIENT ID`` and ``CLIENT SECRET`` into your ``Google Social Account`` on the Admin site.

## Adding Persona Authentication:

The persona authentication is unique as it does not require ids or settings.

* Ensure that ``allauth.socialaccount.providers.persona`` is listed in your ``INSTALLED_APPS``


