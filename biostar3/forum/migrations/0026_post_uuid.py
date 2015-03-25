# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0025_blog_usergroup'),
    ]

    operations = [
        migrations.AddField(
            model_name='post',
            name='uuid',
            field=models.CharField(max_length=256, null=True),
            preserve_default=True,
        ),
    ]
