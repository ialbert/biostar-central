# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0015_field_rename'),
    ]

    operations = [
        migrations.RenameField(
            model_name='groupsub',
            old_name='pref',
            new_name='type',
        ),
        migrations.RenameField(
            model_name='post',
            old_name='group',
            new_name='usergroup',
        ),
        migrations.RenameField(
            model_name='postsub',
            old_name='pref',
            new_name='type',
        ),
        migrations.AlterField(
            model_name='groupperm',
            name='role',
            field=models.IntegerField(default=0, choices=[(0, 'Moderator'), (1, 'Admin')]),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='groupperm',
            unique_together=set([('user', 'usergroup')]),
        ),
        migrations.AlterUniqueTogether(
            name='postsub',
            unique_together=set([('user', 'post')]),
        ),
    ]
