# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import biostar3.forum.models


class Migration(migrations.Migration):

    dependencies = [
        ('taggit', '0001_initial'),
        ('forum', '0005_post_group'),
    ]

    operations = [
        migrations.AddField(
            model_name='post',
            name='tags',
            field=biostar3.forum.models.MyTaggableManager(to='taggit.Tag', through='taggit.TaggedItem', help_text='A comma-separated list of tags.', verbose_name='Tags'),
            preserve_default=True,
        ),
    ]
