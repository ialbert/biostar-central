#
# pyinfra deploy.py file
#
from pyinfra.operations import server, files, apt, git, pip, server

# Destination directory.
DEST = '/tmp/export/www/biostar-central'

# The DJANGO settings module.
DJANGO_SETTINGS_MODULE = 'biostar.forum.settings'

#'conf.run.site_settings'

# Pull from the repository.
git.repo(
    name="Pull from the repository",
    src="git@github.com:ialbert/biostar-central.git",
    dest=f"{DEST}",
    pull=True,
)

files.template(
    name="Create initialization script",
    src="templates/init.sh",
    dest=f"{DEST}/init.sh",
    mode="755",
    DJANGO_SETTINGS_MODULE=f'{DJANGO_SETTINGS_MODULE}',
)

server.shell(
    name="Run initialization script",
    commands=f'cd {DEST} && bash init.sh'
)
