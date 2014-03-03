Install
=======

The sourcecode can be obtained via::

	git clone https://github.com/ialbert/biostar-central.git

Getting started
---------------

Get the source and switch to the source directory. The
recommended installation is via ``virtualenv`` and ``pip``::

	# Install the requirements.
    pip install --upgrade -r conf/requirements/base.txt

	# Load the environment variables.
    source conf/defaults.env

	# Initialize, import test data and run the site.
    ./biostar.sh init import run

Visit ``http://locahost:8080`` to see the site loaded with default settings.

The default admin is ``foo@bar.com`` password ``foobar``. The default email
handler will print to the console. You can reset the password
for any user then copy paste the password reset url into the browser.

Run the manager on its own to see all the commands at your disposal::

	./biostar.sh

To enable searching you must the content with::

    ./biostar.sh index

Social authentication
---------------------

The social logins settings will need to be initialized with the proper
authentication parameters. Typically this involves creating an
application at the provider and obtaining the credentials.

See the ``conf/defaults.env`` for the proper variable naming.

Adding Facebook authentication:

* Create Authentication App: http://developers.facebook.com/setup/
* More information: Facebook Developer Resources: http://developers.facebook.com/docs/authentication/

Adding Google authentication:

* Google Developer Console: https://cloud.google.com/console/project
* Create new project and copy data from credentials
* Callback must be ``http://domain/accounts/google/login/callback/``

Twitter:

* Add your application at Twitter Apps Interface: http://twitter.com/apps/

External authentication
-----------------------

Other domains can provide authentication for Biostar by setting a cookie
with a certain value. For this to work Biostar will have to be set to
run as a subdomain of the hosting site.

Cookie settings
^^^^^^^^^^^^^^^

The cookie value needs to contain the ``email:hash`` as value.
For exampl if the ``EXTERNAL_AUTH`` django settings are::

    # Cookie name, cookie secret key pair
    EXTERNAL_AUTH = [
        ("foo.bar.com", "ABC"),
    ]

If an unauthenticated user sends a cookie named ``foo.bar.com`` with the value::

    foo@bar.com:d46d8c07777e3adf739cfc0c432759b0

then Biostar will automatically log in the user. It will automatically create
an account for the user if the email does not already exist.

Setting the  ``EXTERNAL_LOGIN_URL`` and ``EXTERNAL_LOGOUT_URL`` settings  will also
perform the redirects to the external site login and logout urls::

    EXTERNAL_LOGIN_URL = "http://some.site.com/login"
    EXTERNAL_LOGOUT_URL = "http://some.site.com/logout"

Generating the value is simple like so::

    email = "foo@bar.com"
    digest = hmac.new(key, email).hexdigest()
    value = "%s:%s" % (email, digest)

Prefill post
^^^^^^^^^^^^

Set the ``title``, ``tag_val``, ``content`` and ``category`` fields of a
get request to pre-populate a question::

    http://localhost:8080/p/new/post/?title=Need+help+with+bwa&tag_val=bwa+samtools&content=What+does+it+do?&category=SNP-Calling

Migrating from Biostar 1.X
--------------------------

Due to the complete rework there is no database schema migration.

Instead users of
Biostar 1 site are expected to export their data with a script provided in Biostar 1
then import it with a management command provided with Biostar 2.

The migration will take the following steps:

1. Set the ``BIOSTAR_MIGRATE_DIR`` environment variable to point to a work directory that
   will hold the temporary data, for example  ``export BIOSTAR_MIGRATE_DIR="~/tmp/biostar_export"``

2. Load the environment variables for the Biostar 1 site
   then run ``python -m main.bin.export -u -p -v``. This will dump the contents of the site
   into the directory that ``BIOSTAR_MIGRATE_DIR`` points to.

3. Load the environment variables for you Biostar 2 site then run the
   ``./biostar.sh import_biostar1`` command.

Some caveats, depending how you set the variables you may need to be located in
the root of your site. This applies for the default settings that both sites come
with, as the root is determined relative to the directory that the command is run in.