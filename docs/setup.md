# Site setup and customization

&bull; back to [Documentation Home](index.md)

**Important**: You should not need to directly modify the settings or template files
that come with the source code distribution. Doing so will only make your life
harder later when you need to pull in changes.

The correct approach is to create copies of the files under new names and
use those. We'll run through an example of creating a
new custom site called ``Simple Q&A``.

Please also note that site customization is a different topic than that of [site deployment](deploy.md).

## The settings file

There are example settings files in the ``run`` folder for ``sqlite``, ``postgres`` databases as well
as other settings.

1. Create a directory that will host your site. For example ``mysite``
1. Create a directory that will contain your theme. For exmple ``mysite/theme``
1. Copy ``run/sqlite.env`` and ``run/sqlite.py`` under the names ``mysite/simple.env`` and ``mysite/simple.py``.
1. Create an empty ``__init__.py`` file in the ``mysite`` directory to indicate that it will be used as a package.
1. Edit ``mysite/simple.env`` and set ``DJANGO_SETTINGS_MODULE=mysite.simple``. Change all other settings as necessary.
1. Edit ``mysite/simple.py`` and edit settings that need to be different.

Note that the ``simple.py`` first loads all settings from ``biostar3/settings/base.py`` and
then overrides some settings. For a full list of what can be customized see ``biostar3/settings/base.py``

## The theme files

The customized ``mysite/simple.env`` file should be exporting the variable ``THEME_PATH=mysite/theme``.
For every template Biostar will first look in this directory and if it does not find a template
there falls back to the default location of the source directory. Hence only the templates
that actually change need to be present in the ``mysite/theme`` directory.
