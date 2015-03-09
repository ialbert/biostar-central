# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0010_replytoken'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='usergroup',
            name='author',
        ),
        migrations.AddField(
            model_name='usergroup',
            name='description',
            field=models.TextField(default='default group'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='usergroup',
            name='domain',
            field=models.CharField(default='www', unique=True, max_length=15, db_index=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='usergroup',
            name='owner',
            field=models.ForeignKey(related_name='owners', to=settings.AUTH_USER_MODEL, null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='message_prefs',
            field=models.IntegerField(default=3, choices=[(3, b'Smart mode'), (0, b'Local messages'), (1, b'Email messages'), (4, b'Mail List mode'), (2, b'No Messages')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='usergroup',
            name='name',
            field=models.CharField(unique=True, max_length=25, db_index=True),
            preserve_default=True,
        ),
    ]
