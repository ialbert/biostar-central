# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import biostar3.forum.models


class Migration(migrations.Migration):

    dependencies = [
        ('taggit', '0001_initial'),
        ('forum', '0002_auto_20150625_1755'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='profile',
            name='shortcuts_json',
        ),
        migrations.RemoveField(
            model_name='profile',
            name='shortcuts_text',
        ),
        migrations.AddField(
            model_name='profile',
            name='tags',
            field=biostar3.forum.models.MyTaggableManager(help_text='A comma-separated list of tags.', verbose_name='Tags', through='taggit.TaggedItem', to='taggit.Tag'),
            preserve_default=True,
        ),
    ]
