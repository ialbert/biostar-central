# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0008_digestprefs'),
    ]

    operations = [
        migrations.CreateModel(
            name='GroupPerm',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('type', models.IntegerField(default=0, choices=[(0, 'Write'), (1, 'Admin')])),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='UserGroup',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=15)),
                ('public', models.BooleanField(default=True)),
                ('visible', models.BooleanField(default=True)),
                ('creation_date', models.DateTimeField(auto_now_add=True)),
                ('author', models.ForeignKey(related_name='author', to=settings.AUTH_USER_MODEL, null=True)),
                ('users', models.ManyToManyField(related_name='usergroups', to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.RemoveField(
            model_name='groupinfo',
            name='author',
        ),
        migrations.RemoveField(
            model_name='groupinfo',
            name='group',
        ),
        migrations.DeleteModel(
            name='GroupInfo',
        ),
        migrations.AddField(
            model_name='groupperm',
            name='group',
            field=models.ForeignKey(to='forum.UserGroup'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='groupperm',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AlterModelOptions(
            name='user',
            options={},
        ),
        migrations.AlterField(
            model_name='post',
            name='group',
            field=models.ForeignKey(blank=True, to='forum.UserGroup', null=True),
            preserve_default=True,
        ),
    ]
