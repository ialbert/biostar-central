# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0024_flatpage'),
    ]

    operations = [
        migrations.AddField(
            model_name='blog',
            name='usergroup',
            field=models.ForeignKey(to='forum.UserGroup', null=True),
            preserve_default=True,
        ),
    ]
