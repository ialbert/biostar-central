# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0006_post_tags'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='post',
            name='tag_set',
        ),
        migrations.RemoveField(
            model_name='profile',
            name='tags',
        ),
        migrations.DeleteModel(
            name='Tag',
        ),
    ]
