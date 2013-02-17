External Authentication
=======================

Trusted third party websites may provide user authentication for BioStar. This means
that a user logged into such a website may use their account to post in BioStar under
the same name and have their accounts associated with the originating website.


Conditions/Restrictions
-----------------------

1. A user can only authenticate via the site that they originally registered with.
   In practice this means that the user must log in via the site they established the account through.

2. Authentication may not be transferred. In practice this means that a user
   that first registers via BioStar cannot transfer their registration information to another site and vice versa,
   a user that registers via a third party website may not transfer their
   registration information to BioStar.

3. A user registering via a third party website may not have the same email as an existing BioStar user. 


Trusted Websites
----------------

External authentication takes place via requesting a BioStar URL with an authenticated message in the
URL parameters. The message must be a base64 encoded JSON string accompanied by the HMAC digest created
with a secret key agreed upon by both BioStar and the third party website.

Three parameters need to be passed: `name`, `data` and `digest` corresponding
to the pre-agreed key name, the base64 encoded data and the HMAC digest of the data.

Example
-------

Assuming that the secret key is `abcd`

    name = 'test-key'
    json = { "email":"john.doe@gmail.com", id=1 }
    data = eyAiZW1haWwiOiJqb2huLmRvZUBnbWFpbC5jb20iLCBpZD0xIH0=
    digest =



To log in the user the
