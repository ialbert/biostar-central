# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0003_auto_20150626_1533'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='user',
            name='is_admin',
        ),
        migrations.AlterField(
            model_name='postsub',
            name='type',
            field=models.IntegerField(choices=[(4, 'Smart Mode'), (1, 'Local Messages'), (2, 'Email Messages'), (5, 'Mailing List'), (3, 'No Messages')], default=4),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='message_prefs',
            field=models.IntegerField(choices=[(4, 'Smart Mode'), (1, 'Local Messages'), (2, 'Email Messages'), (5, 'Mailing List'), (3, 'No Messages')], default=4),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='subs_type',
            field=models.IntegerField(choices=[(4, 'Smart Mode'), (1, 'Local Messages'), (2, 'Email Messages'), (5, 'Mailing List'), (3, 'No Messages')], default=4),
            preserve_default=True,
        ),
    ]
