Installation
------------

Get the source and switch to the source directory

The recommended installation is via `virtualenv` and `pip`:

    pip install -r requirements/base.txt

The site manager is `biostar.sh`. This command can take one or more commands like so:

    ./biostar.sh delete init import run

Visit `http://locahost:8080` to see the site loaded with default settings.
The default admin is `foo@bar.com` password `foobar`.
The default email handler will print to the console. You can reset the password
for any user then copy paste the password reset url into the browser.

The Social Logins will need to be enabled via the proper authentication parameters (see `defaults.env`)

To enable searching you must the content with:

    ./biostar.sh index

The next step is to [deploy Biostar][deploy].

[deploy]: docs/deploy.md
