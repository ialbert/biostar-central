# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
from django.utils.timezone import utc


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0020_vote_unread'),
    ]

    operations = [
        migrations.AddField(
            model_name='post',
            name='file',
            field=models.FileField(null=True, upload_to='files/%Y/%m/%d', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='post',
            name='last_activity',
            field=models.DateTimeField(default=datetime.datetime(2015, 3, 21, 13, 7, 13, 965720, tzinfo=utc), db_index=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='user',
            name='portrait',
            field=models.FileField(null=True, upload_to='users/img/%Y/%m/%d', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='usergroup',
            name='css_file',
            field=models.FileField(null=True, upload_to='groups/css/%Y/%m/%d', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='usergroup',
            name='last_activity',
            field=models.DateTimeField(default=datetime.datetime(2015, 3, 21, 13, 7, 27, 851961, tzinfo=utc)),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='usergroup',
            name='post_count',
            field=models.IntegerField(default=0),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='usergroup',
            name='user_count',
            field=models.IntegerField(default=1),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='usergroup',
            name='creation_date',
            field=models.DateTimeField(),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='usergroup',
            name='logo',
            field=models.FileField(null=True, upload_to='groups/img/%Y/%m/%d', blank=True),
            preserve_default=True,
        ),
    ]
