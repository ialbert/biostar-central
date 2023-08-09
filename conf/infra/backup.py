#
# pyinfra backup.py file
#
from pyinfra.operations import server, files, apt, git


backup_src = '/export/www/biostar-central/export/backup/foo.gz'

backup_dest = '/export/www/biostar-central/export/backup/foo.gz'


# Copy over the latest backup file
files.rsync(
    src=backup_src,
    dest=dest,
    flags=['-a',]
)


