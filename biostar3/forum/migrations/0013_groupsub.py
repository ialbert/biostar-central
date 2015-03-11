# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0012_usergroup_logo'),
    ]

    operations = [
        migrations.CreateModel(
            name='GroupSub',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('pref', models.IntegerField(default=3, choices=[(3, b'Default mode'), (1, b'Email mode'), (4, b'Mailing list'), (2, b'No Emails'), (5, b'Leave Group')])),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('usergroup', models.ForeignKey(to='forum.UserGroup')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PostSub',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('pref', models.IntegerField(default=3, choices=[(3, b'Default mode'), (1, b'Email mode'), (4, b'Mailing list'), (2, b'No Emails'), (5, b'Leave Group')])),
                ('post', models.ForeignKey(to='forum.Post')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AlterUniqueTogether(
            name='groupsub',
            unique_together=set([('user', 'usergroup')]),
        ),
        migrations.RemoveField(
            model_name='usergroup',
            name='users',
        ),
        migrations.AlterField(
            model_name='profile',
            name='message_prefs',
            field=models.IntegerField(default=3, choices=[(3, b'Default mode'), (1, b'Email mode'), (4, b'Mailing list'), (2, b'No Emails'), (5, b'Leave Group')]),
            preserve_default=True,
        ),
    ]
