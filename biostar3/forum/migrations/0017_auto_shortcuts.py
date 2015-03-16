# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0016_more_field_renames'),
    ]

    operations = [
        migrations.AddField(
            model_name='profile',
            name='html',
            field=models.TextField(default='', null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='shortcuts_json',
            field=models.TextField(default='', null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='shortcuts_text',
            field=models.TextField(default='', max_length=1000, null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='info',
            field=models.TextField(default='', max_length=5000, null=True, blank=True),
            preserve_default=True,
        ),
    ]
