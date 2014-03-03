# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Tag'
        db.create_table(u'posts_tag', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.TextField')(max_length=50, db_index=True)),
            ('count', self.gf('django.db.models.fields.IntegerField')(default=0)),
        ))
        db.send_create_signal(u'posts', ['Tag'])

        # Adding model 'Post'
        db.create_table(u'posts_post', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('title', self.gf('django.db.models.fields.CharField')(max_length=140)),
            ('author', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['users.User'])),
            ('lastedit_user', self.gf('django.db.models.fields.related.ForeignKey')(related_name=u'editor', to=orm['users.User'])),
            ('rank', self.gf('django.db.models.fields.FloatField')(default=0, blank=True)),
            ('status', self.gf('django.db.models.fields.IntegerField')(default=1)),
            ('type', self.gf('django.db.models.fields.IntegerField')(db_index=True)),
            ('vote_count', self.gf('django.db.models.fields.IntegerField')(default=0, db_index=True, blank=True)),
            ('view_count', self.gf('django.db.models.fields.IntegerField')(default=0, blank=True)),
            ('reply_count', self.gf('django.db.models.fields.IntegerField')(default=0, blank=True)),
            ('comment_count', self.gf('django.db.models.fields.IntegerField')(default=0, blank=True)),
            ('book_count', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('changed', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('subs_count', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('thread_score', self.gf('django.db.models.fields.IntegerField')(default=0, db_index=True, blank=True)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(db_index=True)),
            ('lastedit_date', self.gf('django.db.models.fields.DateTimeField')(db_index=True)),
            ('sticky', self.gf('django.db.models.fields.BooleanField')(default=False, db_index=True)),
            ('has_accepted', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('root', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'descendants', null=True, to=orm['posts.Post'])),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'children', null=True, to=orm['posts.Post'])),
            ('content', self.gf('django.db.models.fields.TextField')(default=u'')),
            ('tag_val', self.gf('django.db.models.fields.CharField')(default=u'', max_length=100, blank=True)),
            ('site', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['sites.Site'], null=True)),
        ))
        db.send_create_signal(u'posts', ['Post'])

        # Adding M2M table for field tag_set on 'Post'
        m2m_table_name = db.shorten_name(u'posts_post_tag_set')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('post', models.ForeignKey(orm[u'posts.post'], null=False)),
            ('tag', models.ForeignKey(orm[u'posts.tag'], null=False))
        ))
        db.create_unique(m2m_table_name, ['post_id', 'tag_id'])

        # Adding model 'PostView'
        db.create_table(u'posts_postview', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('ip', self.gf('django.db.models.fields.GenericIPAddressField')(default=u'', max_length=39, null=True, blank=True)),
            ('post', self.gf('django.db.models.fields.related.ForeignKey')(related_name=u'post_views', to=orm['posts.Post'])),
            ('date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'posts', ['PostView'])

        # Adding model 'Vote'
        db.create_table(u'posts_vote', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('author', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['users.User'])),
            ('post', self.gf('django.db.models.fields.related.ForeignKey')(related_name=u'votes', to=orm['posts.Post'])),
            ('type', self.gf('django.db.models.fields.IntegerField')(db_index=True)),
            ('date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, db_index=True, blank=True)),
        ))
        db.send_create_signal(u'posts', ['Vote'])

        # Adding model 'Subscription'
        db.create_table(u'posts_subscription', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['users.User'])),
            ('post', self.gf('django.db.models.fields.related.ForeignKey')(related_name=u'subs', to=orm['posts.Post'])),
            ('type', self.gf('django.db.models.fields.IntegerField')(default=0, db_index=True)),
            ('date', self.gf('django.db.models.fields.DateTimeField')(db_index=True)),
        ))
        db.send_create_signal(u'posts', ['Subscription'])

        # Adding unique constraint on 'Subscription', fields ['user', 'post']
        db.create_unique(u'posts_subscription', ['user_id', 'post_id'])


    def backwards(self, orm):
        # Removing unique constraint on 'Subscription', fields ['user', 'post']
        db.delete_unique(u'posts_subscription', ['user_id', 'post_id'])

        # Deleting model 'Tag'
        db.delete_table(u'posts_tag')

        # Deleting model 'Post'
        db.delete_table(u'posts_post')

        # Removing M2M table for field tag_set on 'Post'
        db.delete_table(db.shorten_name(u'posts_post_tag_set'))

        # Deleting model 'PostView'
        db.delete_table(u'posts_postview')

        # Deleting model 'Vote'
        db.delete_table(u'posts_vote')

        # Deleting model 'Subscription'
        db.delete_table(u'posts_subscription')


    models = {
        u'posts.post': {
            'Meta': {'object_name': 'Post'},
            'author': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['users.User']"}),
            'book_count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'changed': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'comment_count': ('django.db.models.fields.IntegerField', [], {'default': '0', 'blank': 'True'}),
            'content': ('django.db.models.fields.TextField', [], {'default': "u''"}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'db_index': 'True'}),
            'has_accepted': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'lastedit_date': ('django.db.models.fields.DateTimeField', [], {'db_index': 'True'}),
            'lastedit_user': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'editor'", 'to': u"orm['users.User']"}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'children'", 'null': 'True', 'to': u"orm['posts.Post']"}),
            'rank': ('django.db.models.fields.FloatField', [], {'default': '0', 'blank': 'True'}),
            'reply_count': ('django.db.models.fields.IntegerField', [], {'default': '0', 'blank': 'True'}),
            'root': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'descendants'", 'null': 'True', 'to': u"orm['posts.Post']"}),
            'site': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['sites.Site']", 'null': 'True'}),
            'status': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'sticky': ('django.db.models.fields.BooleanField', [], {'default': 'False', 'db_index': 'True'}),
            'subs_count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'tag_set': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['posts.Tag']", 'symmetrical': 'False', 'blank': 'True'}),
            'tag_val': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '100', 'blank': 'True'}),
            'thread_score': ('django.db.models.fields.IntegerField', [], {'default': '0', 'db_index': 'True', 'blank': 'True'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '140'}),
            'type': ('django.db.models.fields.IntegerField', [], {'db_index': 'True'}),
            'view_count': ('django.db.models.fields.IntegerField', [], {'default': '0', 'blank': 'True'}),
            'vote_count': ('django.db.models.fields.IntegerField', [], {'default': '0', 'db_index': 'True', 'blank': 'True'})
        },
        u'posts.postview': {
            'Meta': {'object_name': 'PostView'},
            'date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'ip': ('django.db.models.fields.GenericIPAddressField', [], {'default': "u''", 'max_length': '39', 'null': 'True', 'blank': 'True'}),
            'post': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'post_views'", 'to': u"orm['posts.Post']"})
        },
        u'posts.subscription': {
            'Meta': {'unique_together': "((u'user', u'post'),)", 'object_name': 'Subscription'},
            'date': ('django.db.models.fields.DateTimeField', [], {'db_index': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'post': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'subs'", 'to': u"orm['posts.Post']"}),
            'type': ('django.db.models.fields.IntegerField', [], {'default': '0', 'db_index': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['users.User']"})
        },
        u'posts.tag': {
            'Meta': {'object_name': 'Tag'},
            'count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {'max_length': '50', 'db_index': 'True'})
        },
        u'posts.vote': {
            'Meta': {'object_name': 'Vote'},
            'author': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['users.User']"}),
            'date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'db_index': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'post': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'votes'", 'to': u"orm['posts.Post']"}),
            'type': ('django.db.models.fields.IntegerField', [], {'db_index': 'True'})
        },
        u'sites.site': {
            'Meta': {'ordering': "(u'domain',)", 'object_name': 'Site', 'db_table': "u'django_site'"},
            'domain': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'users.user': {
            'Meta': {'object_name': 'User'},
            'badges': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'email': ('django.db.models.fields.EmailField', [], {'unique': 'True', 'max_length': '255', 'db_index': 'True'}),
            'flair': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '15'}),
            'full_score': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_admin': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'name': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '255'}),
            'new_messages': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'score': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'site': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['sites.Site']", 'null': 'True'}),
            'status': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'type': ('django.db.models.fields.IntegerField', [], {'default': '0'})
        }
    }

    complete_apps = ['posts']