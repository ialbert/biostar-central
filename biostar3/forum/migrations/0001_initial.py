# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='User',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('password', models.CharField(max_length=128, verbose_name='password')),
                ('last_login', models.DateTimeField(default=django.utils.timezone.now, verbose_name='last login')),
                ('email', models.EmailField(unique=True, max_length=255, verbose_name=b'Email', db_index=True)),
                ('name', models.CharField(default=b'Biostar User', max_length=255, verbose_name=b'Name')),
                ('is_active', models.BooleanField(default=True)),
                ('is_admin', models.BooleanField(default=False)),
                ('is_staff', models.BooleanField(default=False)),
                ('type', models.IntegerField(default=0, choices=[(0, b'User'), (1, b'Moderator'), (2, b'Admin'), (3, b'Blog')])),
                ('status', models.IntegerField(default=0, choices=[(0, b'New User'), (1, b'Trusted'), (2, b'Suspended'), (3, b'Banned')])),
                ('new_messages', models.IntegerField(default=0)),
                ('badges', models.IntegerField(default=0)),
                ('score', models.IntegerField(default=0)),
                ('activity', models.IntegerField(default=0)),
                ('flair', models.CharField(default=b'', max_length=15, verbose_name=b'Flair')),
            ],
            options={
                'db_table': 'users_user',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Post',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('title', models.CharField(max_length=200)),
                ('rank', models.FloatField(default=0, blank=True)),
                ('status', models.IntegerField(default=1, choices=[(0, b'Pending'), (1, b'Open'), (2, b'Closed'), (3, b'Deleted')])),
                ('type', models.IntegerField(db_index=True, choices=[(0, b'Question'), (1, b'Answer'), (6, b'Comment'), (2, b'Job'), (3, b'Forum'), (8, b'Tutorial'), (7, b'Data'), (4, b'Page'), (10, b'Tool'), (11, b'News'), (5, b'Blog'), (9, b'Bulletin Board')])),
                ('vote_count', models.IntegerField(default=0, db_index=True, blank=True)),
                ('view_count', models.IntegerField(default=0, blank=True)),
                ('reply_count', models.IntegerField(default=0, blank=True)),
                ('comment_count', models.IntegerField(default=0, blank=True)),
                ('book_count', models.IntegerField(default=0)),
                ('changed', models.BooleanField(default=True)),
                ('subs_count', models.IntegerField(default=0)),
                ('thread_score', models.IntegerField(default=0, db_index=True, blank=True)),
                ('creation_date', models.DateTimeField(db_index=True)),
                ('lastedit_date', models.DateTimeField(db_index=True)),
                ('sticky', models.BooleanField(default=False, db_index=True)),
                ('has_accepted', models.BooleanField(default=False)),
                ('content', models.TextField(default=b'')),
                ('html', models.TextField(default=b'')),
                ('tag_val', models.CharField(default=b'', max_length=100, blank=True)),
                ('author', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('lastedit_user', models.ForeignKey(related_name='editor', to=settings.AUTH_USER_MODEL)),
                ('parent', models.ForeignKey(related_name='children', blank=True, to='forum.Post', null=True)),
                ('root', models.ForeignKey(related_name='descendants', blank=True, to='forum.Post', null=True)),
            ],
            options={
                'db_table': 'posts_post',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Tag',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.TextField(max_length=50, db_index=True)),
                ('count', models.IntegerField(default=0)),
            ],
            options={
                'db_table': 'posts_tag',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='post',
            name='tag_set',
            field=models.ManyToManyField(to='forum.Tag', blank=True),
            preserve_default=True,
        ),
    ]
