# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0007_delete_old_tag'),
    ]

    operations = [
        migrations.AddField(
            model_name='groupinfo',
            name='public',
            field=models.BooleanField(default=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='digest_prefs',
            field=models.IntegerField(default=2, choices=[(0, 'Never'), (1, 'Daily'), (2, 'Weekly'), (3, 'Monthly')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='post',
            name='lastedit_user',
            field=models.ForeignKey(related_name='editor', to=settings.AUTH_USER_MODEL, null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='flair',
            field=models.CharField(default='', max_length=15, verbose_name='Flair', blank=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='name',
            field=models.CharField(default='', max_length=255, verbose_name='Name'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='site',
            field=models.ForeignKey(blank=True, to='sites.Site', null=True),
            preserve_default=True,
        ),
    ]
