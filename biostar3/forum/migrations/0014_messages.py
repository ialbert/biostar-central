# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('forum', '0013_groupsub'),
    ]

    operations = [
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
        migrations.AlterField(
            model_name='groupsub',
            name='pref',
            field=models.IntegerField(default=3, choices=[(3, b'Smart mode'), (2, b'Local messages'), (0, b'Local tracker'), (1, b'Email tracker'), (4, b'Mailing list'), (5, b'Leave Group')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='postsub',
            name='pref',
            field=models.IntegerField(default=3, choices=[(3, b'Smart mode'), (2, b'Local messages'), (0, b'Local tracker'), (1, b'Email tracker'), (4, b'Mailing list'), (5, b'Leave Group')]),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='profile',
            name='message_prefs',
            field=models.IntegerField(default=3, choices=[(3, b'Smart mode'), (2, b'Local messages'), (0, b'Local tracker'), (1, b'Email tracker'), (4, b'Mailing list'), (5, b'Leave Group')]),
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
