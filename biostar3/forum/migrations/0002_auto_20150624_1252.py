# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.utils.timezone import utc
import datetime
from django.conf import settings
import biostar3.forum.models


class Migration(migrations.Migration):

    dependencies = [
        ('taggit', '0001_initial'),
        ('auth', '0001_initial'),
        ('forum', '0001_initial'),
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
                ('uuid', models.CharField(default='', unique=True, max_length=100, blank=True)),
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
        migrations.CreateModel(
            name='FederatedContent',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('obj_id', models.IntegerField(default=0, db_index=True)),
                ('domain', models.TextField(default='')),
                ('content', models.TextField(default='')),
                ('changed', models.BooleanField(default=False)),
                ('creation_date', models.DateTimeField(auto_now=True, db_index=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='FlatPage',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
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
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('unread', models.BooleanField(default=True)),
                ('date', models.DateTimeField(null=True, db_index=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='MessageBody',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
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
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('type', models.IntegerField(default=0, choices=[(0, b'Smart Mode'), (1, b'Local Tracker'), (2, b'Email Tracker'), (3, b'Mailing List'), (4, b'Leave Group')])),
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
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
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
            name='TaggedUser',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
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
            field=models.ForeignKey(blank=True, to='forum.Post', null=True),
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
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 24, 12, 52, 36, 220302, tzinfo=utc), db_index=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='post',
            name='tags',
            field=biostar3.forum.models.MyTaggableManager(to='taggit.Tag', through='taggit.TaggedItem', help_text='A comma-separated list of tags.', verbose_name='Tags'),
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
            field=models.IntegerField(default=2, choices=[(0, b'Never'), (1, b'Daily'), (2, b'Weekly'), (3, b'Monthly')]),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='html',
            field=models.TextField(default='', null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='shortcuts_json',
            field=models.TextField(default='', null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='profile',
            name='shortcuts_text',
            field=models.TextField(default='', max_length=1000, null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='groups',
            field=models.ManyToManyField(related_query_name='user', related_name='user_set', to='auth.Group', blank=True, help_text='The groups this user belongs to. A user will get all permissions granted to each of his/her group.', verbose_name='groups'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='handle',
            field=models.CharField(default='', max_length=25, verbose_name='Handle', blank=True),
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
            field=models.FileField(null=True, upload_to='users/img/%Y/%m/%d', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='subs_type',
            field=models.IntegerField(default=0, choices=[(0, b'Smart Mode'), (1, b'Local Tracker'), (2, b'Email Tracker'), (3, b'Mailing List'), (4, b'Leave Group')]),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='user',
            name='user_permissions',
            field=models.ManyToManyField(related_query_name='user', related_name='user_set', to='auth.Permission', blank=True, help_text='Specific permissions for this user.', verbose_name='user permissions'),
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
            field=models.ForeignKey(related_name='editor', to=settings.AUTH_USER_MODEL, null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='info',
            field=models.TextField(default='', max_length=5000, null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='message_prefs',
            field=models.IntegerField(default=0, choices=[(0, b'Smart Mode'), (1, b'Local Tracker'), (2, b'Email Tracker'), (3, b'Mailing List'), (4, b'Leave Group')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='flair',
            field=models.CharField(default='0,0,0,0', max_length=255, verbose_name='Flair', blank=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='name',
            field=models.CharField(default='', max_length=255, verbose_name='Name'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='user',
            name='site',
            field=models.ForeignKey(blank=True, to='sites.Site', null=True),
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
