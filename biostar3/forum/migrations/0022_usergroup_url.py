# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0021_more_filefields'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='usergroup',
            name='description',
        ),
        migrations.AddField(
            model_name='usergroup',
            name='html',
            field=models.TextField(default='group info turned into html'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='usergroup',
            name='info',
            field=models.TextField(default='default group info'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='usergroup',
            name='url',
            field=models.CharField(default='full url to the group home', max_length=255),
            preserve_default=True,
        ),
    ]
