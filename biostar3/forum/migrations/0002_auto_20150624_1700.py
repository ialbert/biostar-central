# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.utils.timezone import utc
from django.conf import settings
import biostar3.forum.models
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('taggit', '0001_initial'),
        ('sites', '0001_initial'),
        ('auth', '0001_initial'),
        ('forum', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Award',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('date', models.DateTimeField()),
                ('context', models.CharField(max_length=1000, default='')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Badge',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('uuid', models.CharField(max_length=100, unique=True, blank=True, default='')),
                ('name', models.CharField(max_length=50)),
                ('desc', models.CharField(max_length=200, default='')),
                ('type', models.IntegerField(default=0, choices=[(0, 'User badge'), (1, 'Post badge')])),
                ('style', models.IntegerField(default=0, choices=[(0, 'Bronze'), (1, 'Silver'), (2, 'Gold')])),
                ('unique', models.BooleanField(default=False)),
                ('count', models.IntegerField(default=0)),
                ('icon', models.CharField(max_length=250, default='<i class="fa fa-asterisk"></i>')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='FederatedContent',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('obj_id', models.IntegerField(db_index=True, default=0)),
                ('domain', models.TextField(default='')),
                ('content', models.TextField(default='')),
                ('changed', models.BooleanField(default=False)),
                ('creation_date', models.DateTimeField(db_index=True, auto_now=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='FlatPage',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('slug', models.SlugField(default='slug')),
                ('post', models.ForeignKey(to='forum.Post')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Message',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('unread', models.BooleanField(default=True)),
                ('date', models.DateTimeField(db_index=True, null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='MessageBody',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('subject', models.CharField(max_length=500)),
                ('html', models.TextField()),
                ('content', models.TextField()),
                ('date', models.DateTimeField()),
                ('author', models.ForeignKey(related_name='message_bodies', to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PostSub',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('type', models.IntegerField(default=0, choices=[(0, 'Smart Mode'), (1, 'Local Tracker'), (2, 'Email Tracker'), (3, 'Mailing List')])),
                ('post', models.ForeignKey(to='forum.Post')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Post Subscription',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ReplyToken',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('date', models.DateTimeField(auto_created=True)),
                ('token', models.CharField(max_length=256)),
                ('post', models.ForeignKey(to='forum.Post')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SiteSub',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('site', models.ForeignKey(to='sites.Site')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TaggedUser',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('post', models.ForeignKey(to='forum.Post')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AlterUniqueTogether(
            name='postsub',
            unique_together=set([('user', 'post')]),
        ),
        migrations.AddField(
            model_name='message',
            name='body',
            field=models.ForeignKey(related_name='messages', to='forum.MessageBody'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='message',
            name='user',
            field=models.ForeignKey(related_name='messages', to=settings.AUTH_USER_MODEL),
            preserve_default=True,
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
            field=models.ForeignKey(blank=True, null=True, to='forum.Post'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='award',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AlterModelOptions(
            name='post',
            options={'ordering': ['-lastedit_date'], 'permissions': (('moderate_post', 'Can moderate a post'),)},
        ),
        migrations.AlterModelOptions(
            name='postview',
            options={'verbose_name': 'Post View'},
        ),
        migrations.RemoveField(
            model_name='post',
            name='tag_set',
        ),
        migrations.RemoveField(
            model_name='profile',
            name='tags',
        ),
        migrations.DeleteModel(
            name='Tag',
        ),
        migrations.AddField(
            model_name='post',
            name='file',
            field=models.FileField(null=True, upload_to='files/%Y/%m/%d', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='post',
            name='last_activity',
            field=models.DateTimeField(db_index=True, default=datetime.datetime(2015, 6, 24, 17, 0, 48, 517473, tzinfo=utc)),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='post',
            name='tags',
            field=biostar3.forum.models.MyTaggableManager(verbose_name='Tags', help_text='A comma-separated list of tags.', through='taggit.TaggedItem', to='taggit.Tag'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='post',
            name='uuid',
            field=models.CharField(null=True, max_length=256),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='digest_prefs',
            field=models.IntegerField(default=2, choices=[(0, 'Never'), (1, 'Daily'), (2, 'Weekly'), (3, 'Monthly')]),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='html',
            field=models.TextField(null=True, default='', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='shortcuts_json',
            field=models.TextField(null=True, default='', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='shortcuts_text',
            field=models.TextField(null=True, max_length=1000, blank=True, default=''),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='groups',
            field=models.ManyToManyField(related_query_name='user', verbose_name='groups', help_text='The groups this user belongs to. A user will get all permissions granted to each of his/her group.', blank=True, related_name='user_set', to='auth.Group'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='handle',
            field=models.CharField(max_length=25, verbose_name='Handle', blank=True, default=''),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='is_superuser',
            field=models.BooleanField(default=False, verbose_name='superuser status', help_text='Designates that this user has all permissions without explicitly assigning them.'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='portrait',
            field=models.FileField(null=True, upload_to='users/img/%Y/%m/%d', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='subs_type',
            field=models.IntegerField(default=0, choices=[(0, 'Smart Mode'), (1, 'Local Tracker'), (2, 'Email Tracker'), (3, 'Mailing List')]),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='user_permissions',
            field=models.ManyToManyField(related_query_name='user', verbose_name='user permissions', help_text='Specific permissions for this user.', blank=True, related_name='user_set', to='auth.Permission'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='vote',
            name='unread',
            field=models.BooleanField(default=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='post',
            name='lastedit_user',
            field=models.ForeignKey(null=True, related_name='editor', to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='info',
            field=models.TextField(null=True, max_length=5000, blank=True, default=''),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='message_prefs',
            field=models.IntegerField(default=0, choices=[(0, 'Smart Mode'), (1, 'Local Tracker'), (2, 'Email Tracker'), (3, 'Mailing List')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='flair',
            field=models.CharField(max_length=255, verbose_name='Flair', blank=True, default='0,0,0,0'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='name',
            field=models.CharField(max_length=255, verbose_name='Name', default=''),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='site',
            field=models.ForeignKey(blank=True, null=True, to='sites.Site'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='vote',
            name='date',
            field=models.DateTimeField(db_index=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='vote',
            name='type',
            field=models.IntegerField(choices=[(0, 'Upvote'), (1, 'DownVote'), (2, 'Bookmark'), (3, 'Accept')]),
            preserve_default=True,
        ),
    ]
