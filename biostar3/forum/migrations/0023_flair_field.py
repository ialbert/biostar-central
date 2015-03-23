# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0022_usergroup_url'),
    ]

    operations = [
        migrations.AlterField(
            model_name='user',
            name='flair',
            field=models.CharField(default='0,0,0,0', max_length=255, verbose_name='Flair', blank=True),
            preserve_default=True,
        ),
    ]
