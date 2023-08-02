#
# pyinfra deploy.py file
#
from pyinfra.operations import server, files, apt, git

# Pull from the repository.
git.repo(
    name="Pull from the repository",
    src="git@github.com:ialbert/biostar-central.git",
    dest="/export/www/biostar-central",
    pull = True,
)



