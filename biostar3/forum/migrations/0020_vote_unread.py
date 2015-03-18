# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0019_field_changes'),
    ]

    operations = [
        migrations.AddField(
            model_name='vote',
            name='unread',
            field=models.BooleanField(default=True),
            preserve_default=True,
        ),
    ]
