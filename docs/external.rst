
External Authentication
=======================

Other domains can provide authentication for Biostar by setting a cookie
with a certain value. For this to work Biostar will have to be set to
run as a subdomain of the hosting site.

Cookie Settings
---------------

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
------------

Set the ``title``, ``tag_val``, ``content`` and ``category`` fields of a
get request to pre-populate a question::

    http://localhost:8080/p/new/post/?title=Need+help+with+bwa&tag_val=bwa+samtools&content=What+does+it+do?&category=SNP-Calling