# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0011_groups'),
    ]

    operations = [
        migrations.AddField(
            model_name='usergroup',
            name='logo',
            field=models.FileField(null=True, upload_to='groups', blank=True),
            preserve_default=True,
        ),
    ]
