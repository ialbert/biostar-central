This directory is for convenience only.

We put the settings file used in deployment into this directory.

In addition biostar.forum.settings module will also attempts to import

conf/run/secrets.py

module in this directory. We use this module to store to store various sensitive settings: Google Login tokens,
ADMIN passwords etc. that are not public knowledge.

It is not a problem if you don't have that module here.
