# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0017_auto_shortcuts'),
    ]

    operations = [
        migrations.CreateModel(
            name='Award',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('date', models.DateTimeField()),
                ('context', models.CharField(default='', max_length=1000)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Badge',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('desc', models.CharField(default='', max_length=200)),
                ('type', models.IntegerField(default=0, choices=[(0, 'User badge'), (1, 'Post badge')])),
                ('style', models.IntegerField(default=0, choices=[(0, 'Bronze'), (1, 'Silver'), (2, 'Gold')])),
                ('unique', models.BooleanField(default=False)),
                ('count', models.IntegerField(default=0)),
                ('icon', models.CharField(default='<i class="fa fa-asterisk"></i>', max_length=250)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='award',
            name='badge',
            field=models.ForeignKey(to='forum.Badge'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='award',
            name='post',
            field=models.ForeignKey(blank=True, to='forum.Post', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='award',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
    ]
