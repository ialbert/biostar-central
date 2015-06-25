# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.utils.timezone import utc
from django.conf import settings
import biostar3.forum.models
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('sites', '0001_initial'),
        ('auth', '0001_initial'),
        ('taggit', '0001_initial'),
        ('forum', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Award',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
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
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('uuid', models.CharField(blank=True, max_length=100, unique=True, default='')),
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
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('obj_id', models.IntegerField(default=0, db_index=True)),
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
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
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
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
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
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('subject', models.CharField(max_length=500)),
                ('html', models.TextField()),
                ('content', models.TextField()),
                ('date', models.DateTimeField()),
                ('author', models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='message_bodies')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PostSub',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('type', models.IntegerField(default=0, choices=[(0, 'Smart Mode'), (1, 'Local Notifications'), (2, 'Email Notifications'), (3, 'Mailing List')])),
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
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
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
            name='SiteLog',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('date', models.DateTimeField(db_index=True, null=True)),
                ('text', models.TextField()),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL, null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SiteSub',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('site', models.ForeignKey(to='sites.Site')),
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
            field=models.ForeignKey(to='forum.MessageBody', related_name='messages'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='message',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='messages'),
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
            field=models.ForeignKey(null=True, blank=True, to='forum.Post'),
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
            field=models.FileField(blank=True, upload_to='files/%Y/%m/%d', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='post',
            name='last_activity',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 25, 17, 55, 38, 736768, tzinfo=utc), db_index=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='post',
            name='tags',
            field=biostar3.forum.models.MyTaggableManager(to='taggit.Tag', verbose_name='Tags', help_text='A comma-separated list of tags.', through='taggit.TaggedItem'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='post',
            name='uuid',
            field=models.CharField(max_length=256, null=True),
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
            field=models.TextField(blank=True, default='', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='shortcuts_json',
            field=models.TextField(blank=True, default='', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='shortcuts_text',
            field=models.TextField(blank=True, max_length=1000, null=True, default=''),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='groups',
            field=models.ManyToManyField(help_text='The groups this user belongs to. A user will get all permissions granted to each of his/her group.', related_name='user_set', blank=True, related_query_name='user', to='auth.Group', verbose_name='groups'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='handle',
            field=models.CharField(blank=True, max_length=25, verbose_name='Handle', default=''),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='is_superuser',
            field=models.BooleanField(default=False, help_text='Designates that this user has all permissions without explicitly assigning them.', verbose_name='superuser status'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='portrait',
            field=models.FileField(blank=True, upload_to='users/img/%Y/%m/%d', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='subs_type',
            field=models.IntegerField(default=0, choices=[(0, 'Smart Mode'), (1, 'Local Notifications'), (2, 'Email Notifications'), (3, 'Mailing List')]),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='user_permissions',
            field=models.ManyToManyField(help_text='Specific permissions for this user.', related_name='user_set', blank=True, related_query_name='user', to='auth.Permission', verbose_name='user permissions'),
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
            field=models.TextField(blank=True, max_length=5000, null=True, default=''),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='message_prefs',
            field=models.IntegerField(default=0, choices=[(0, 'Smart Mode'), (1, 'Local Notifications'), (2, 'Email Notifications'), (3, 'Mailing List')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='flair',
            field=models.CharField(blank=True, max_length=255, verbose_name='Flair', default='0,0,0,0'),
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
            field=models.ForeignKey(null=True, blank=True, to='sites.Site'),
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
