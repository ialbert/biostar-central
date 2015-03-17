# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('sites', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='User',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('password', models.CharField(max_length=128, verbose_name='password')),
                ('last_login', models.DateTimeField(default=django.utils.timezone.now, verbose_name='last login')),
                ('email', models.EmailField(unique=True, max_length=255, verbose_name='Email', db_index=True)),
                ('name', models.CharField(default='Biostar User', max_length=255, verbose_name='Name')),
                ('is_active', models.BooleanField(default=True)),
                ('is_admin', models.BooleanField(default=False)),
                ('is_staff', models.BooleanField(default=False)),
                ('type', models.IntegerField(default=0, choices=[(0, 'User'), (1, 'Moderator'), (2, 'Admin'), (3, 'Blog')])),
                ('status', models.IntegerField(default=0, choices=[(0, 'New User'), (1, 'Trusted'), (2, 'Suspended'), (3, 'Banned')])),
                ('new_messages', models.IntegerField(default=0)),
                ('badges', models.IntegerField(default=0)),
                ('score', models.IntegerField(default=0)),
                ('activity', models.IntegerField(default=0)),
                ('flair', models.CharField(default='', max_length=15, verbose_name='Flair')),
                ('site', models.ForeignKey(to='sites.Site', null=True)),
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
                ('title', models.CharField(max_length=250)),
                ('rank', models.FloatField(default=0, blank=True)),
                ('status', models.IntegerField(default=1, choices=[(0, 'Pending'), (1, 'Open'), (2, 'Closed'), (3, 'Deleted')])),
                ('type', models.IntegerField(db_index=True, choices=[(0, 'Question'), (1, 'Answer'), (6, 'Comment'), (2, 'Job'), (3, 'Forum'), (8, 'Tutorial'), (7, 'Data'), (4, 'Page'), (10, 'Tool'), (11, 'News'), (5, 'Blog'), (9, 'Bulletin Board')])),
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
                ('content', models.TextField(default='')),
                ('html', models.TextField(default='')),
                ('tag_val', models.CharField(default='', max_length=100, blank=True)),
                ('author', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('lastedit_user', models.ForeignKey(related_name='editor', to=settings.AUTH_USER_MODEL)),
                ('parent', models.ForeignKey(related_name='children', blank=True, to='forum.Post', null=True)),
                ('root', models.ForeignKey(related_name='descendants', blank=True, to='forum.Post', null=True)),
                ('site', models.ForeignKey(to='sites.Site', null=True)),
            ],
            options={
                'ordering': ['-lastedit_date'],
                'db_table': 'posts_post',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PostView',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('ip', models.GenericIPAddressField(default='', null=True, blank=True)),
                ('date', models.DateTimeField(auto_now=True)),
                ('post', models.ForeignKey(related_name='post_views', to='forum.Post')),
            ],
            options={
                'db_table': 'posts_postview',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Profile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('uuid', models.CharField(unique=True, max_length=255, db_index=True)),
                ('last_login', models.DateTimeField()),
                ('date_joined', models.DateTimeField()),
                ('location', models.CharField(default='', max_length=255, blank=True)),
                ('website', models.URLField(default='', max_length=255, blank=True)),
                ('scholar', models.CharField(default='', max_length=255, blank=True)),
                ('twitter_id', models.CharField(default='', max_length=255, blank=True)),
                ('my_tags', models.TextField(default='', max_length=255, blank=True)),
                ('info', models.TextField(default='', null=True, blank=True)),
                ('message_prefs', models.IntegerField(default=3, choices=[(3, b'smart mode'), (0, b'local messages'), (1, b'emails'), (4, b'email for every new thread (mailing list mode)')])),
                ('flag', models.IntegerField(default=0)),
                ('watched_tags', models.CharField(default='', max_length=250, blank=True)),
            ],
            options={
                'db_table': 'users_profile',
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
        migrations.CreateModel(
            name='Vote',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('type', models.IntegerField(db_index=True, choices=[(0, 'Upvote'), (1, 'DownVote'), (2, 'Bookmark'), (3, 'Accept')])),
                ('date', models.DateTimeField(auto_now=True, db_index=True)),
                ('author', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('post', models.ForeignKey(related_name='votes', to='forum.Post')),
            ],
            options={
                'db_table': 'posts_vote',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='profile',
            name='tags',
            field=models.ManyToManyField(to='forum.Tag', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='user',
            field=models.OneToOneField(to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='post',
            name='tag_set',
            field=models.ManyToManyField(to='forum.Tag', blank=True),
            preserve_default=True,
        ),
        migrations.CreateModel(
            name='Blog',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('title', models.CharField(default='', max_length=255, verbose_name='Blog Name')),
                ('desc', models.TextField(default='', blank=True)),
                ('feed', models.URLField()),
                ('link', models.URLField()),
                ('active', models.BooleanField(default=True)),
                ('list_order', models.IntegerField(default=0)),
            ],
            options={
                'db_table': 'planet_blog',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='BlogPost',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('uid', models.CharField(default='', max_length=200)),
                ('title', models.CharField(max_length=200)),
                ('content', models.TextField(default='', max_length=20000)),
                ('html', models.TextField(default='')),
                ('creation_date', models.DateTimeField(db_index=True)),
                ('insert_date', models.DateTimeField(null=True, db_index=True)),
                ('published', models.BooleanField(default=False)),
                ('link', models.URLField()),
                ('blog', models.ForeignKey(to='forum.Blog')),
            ],
            options={
                'db_table': 'planet_blogpost',
            },
            bases=(models.Model,),
        ),
    ]
