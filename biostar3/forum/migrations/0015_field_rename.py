# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0014_messages'),
    ]

    operations = [
        migrations.RenameField(
            model_name='groupperm',
            old_name='type',
            new_name='role',
        ),
        migrations.RenameField(
            model_name='groupperm',
            old_name='group',
            new_name='usergroup',
        ),
    ]
