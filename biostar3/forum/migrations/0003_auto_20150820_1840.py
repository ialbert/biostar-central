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
        migrations.RemoveField(
            model_name='user',
            name='handle',
        ),
        migrations.RemoveField(
            model_name='user',
            name='is_admin',
        ),
        migrations.AddField(
            model_name='message',
            name='post',
            field=models.ForeignKey(to='forum.Post', related_name='messages', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='tags',
            field=biostar3.forum.models.MyTaggableManager(verbose_name='Tags', through='taggit.TaggedItem', to='taggit.Tag', help_text='A comma-separated list of tags.'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='postsub',
            name='type',
            field=models.IntegerField(choices=[(3, 'Smart Mode'), (0, 'Local Messages'), (1, 'Email Messages'), (4, 'Mailing List'), (2, 'No Messages')], default=3),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='message_prefs',
            field=models.IntegerField(choices=[(3, 'Smart Mode'), (0, 'Local Messages'), (1, 'Email Messages'), (4, 'Mailing List'), (2, 'No Messages')], default=3),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='subs_type',
            field=models.IntegerField(choices=[(3, 'Smart Mode'), (0, 'Local Messages'), (1, 'Email Messages'), (4, 'Mailing List'), (2, 'No Messages')], default=3),
            preserve_default=True,
        ),
    ]
