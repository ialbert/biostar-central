# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0002_federatedcontent'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='post',
            options={'ordering': ['-lastedit_date'], 'permissions': (('moderate_post', 'Can moderate a post'),)},
        ),
        migrations.AlterModelOptions(
            name='user',
            options={'permissions': (('moderate_user', 'Can moderate a user'), ('ban_user', 'Can ban a user'))},
        ),
    ]
