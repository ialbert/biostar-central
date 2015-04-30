# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0026_post_uuid'),
    ]

    operations = [
        migrations.AddField(
            model_name='user',
            name='handle',
            field=models.CharField(default='', max_length=25, blank=True, verbose_name='Handle'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='groupsub',
            name='type',
            field=models.IntegerField(default=3, choices=[(3, 'Smart mode'), (2, 'Local messages'), (0, 'Local tracker'), (1, 'Email tracker'), (4, 'Mailing list'), (5, 'Leave Group')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='postsub',
            name='type',
            field=models.IntegerField(default=3, choices=[(3, 'Smart mode'), (2, 'Local messages'), (0, 'Local tracker'), (1, 'Email tracker'), (4, 'Mailing list'), (5, 'Leave Group')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='message_prefs',
            field=models.IntegerField(default=3, choices=[(3, 'Smart mode'), (2, 'Local messages'), (0, 'Local tracker'), (1, 'Email tracker'), (4, 'Mailing list'), (5, 'Leave Group')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='usergroup',
            name='css_file',
            field=models.FileField(upload_to='groups/css', null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='usergroup',
            name='logo',
            field=models.FileField(upload_to='groups/img', null=True, blank=True),
            preserve_default=True,
        ),
    ]
