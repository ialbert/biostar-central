# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'UserProfile'
        db.create_table('server_userprofile', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.OneToOneField')(related_name='profile', unique=True, to=orm['auth.User'])),
            ('display_name', self.gf('django.db.models.fields.CharField')(default='User', max_length=250, db_index=True)),
            ('type', self.gf('django.db.models.fields.IntegerField')(default=1)),
            ('uuid', self.gf('django.db.models.fields.TextField')(unique=True, db_index=True)),
            ('score', self.gf('django.db.models.fields.IntegerField')(default=0, blank=True)),
            ('bronze_badges', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('silver_badges', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('gold_badges', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('new_messages', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('last_visited', self.gf('django.db.models.fields.DateTimeField')()),
            ('status', self.gf('django.db.models.fields.IntegerField')(default=1)),
            ('about_me', self.gf('django.db.models.fields.TextField')(default='', null=True, blank=True)),
            ('about_me_html', self.gf('django.db.models.fields.TextField')(default='', null=True, blank=True)),
            ('location', self.gf('django.db.models.fields.TextField')(default='', null=True, blank=True)),
            ('website', self.gf('django.db.models.fields.URLField')(default='', max_length=250, null=True, blank=True)),
            ('my_tags', self.gf('django.db.models.fields.TextField')(default='', max_length=250, null=True, blank=True)),
            ('scholar', self.gf('django.db.models.fields.TextField')(default='', max_length=50, null=True, blank=True)),
        ))
        db.send_create_signal('server', ['UserProfile'])

        # Adding model 'Tag'
        db.create_table('server_tag', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.TextField')(max_length=50, db_index=True)),
            ('count', self.gf('django.db.models.fields.IntegerField')(default=0)),
        ))
        db.send_create_signal('server', ['Tag'])

        # Adding model 'Post'
        db.create_table('server_post', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('author', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('content', self.gf('django.db.models.fields.TextField')(max_length=10000)),
            ('html', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('title', self.gf('django.db.models.fields.TextField')(max_length=200)),
            ('slug', self.gf('django.db.models.fields.SlugField')(max_length=200, blank=True)),
            ('tag_val', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('views', self.gf('django.db.models.fields.IntegerField')(default=0, db_index=True, blank=True)),
            ('score', self.gf('django.db.models.fields.IntegerField')(default=0, db_index=True, blank=True)),
            ('full_score', self.gf('django.db.models.fields.IntegerField')(default=0, db_index=True, blank=True)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(db_index=True)),
            ('lastedit_date', self.gf('django.db.models.fields.DateTimeField')()),
            ('lastedit_user', self.gf('django.db.models.fields.related.ForeignKey')(related_name='editor', to=orm['auth.User'])),
            ('changed', self.gf('django.db.models.fields.BooleanField')(default=True, db_index=True)),
            ('status', self.gf('django.db.models.fields.IntegerField')(default=100)),
            ('type', self.gf('django.db.models.fields.IntegerField')(db_index=True)),
            ('root', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='descendants', null=True, to=orm['server.Post'])),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='children', null=True, to=orm['server.Post'])),
            ('answer_count', self.gf('django.db.models.fields.IntegerField')(default=0, blank=True)),
            ('accepted', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('url', self.gf('django.db.models.fields.URLField')(default='', max_length=200, blank=True)),
            ('rank', self.gf('django.db.models.fields.FloatField')(default=0, blank=True)),
        ))
        db.send_create_signal('server', ['Post'])

        # Adding M2M table for field tag_set on 'Post'
        db.create_table('server_post_tag_set', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('post', models.ForeignKey(orm['server.post'], null=False)),
            ('tag', models.ForeignKey(orm['server.tag'], null=False))
        ))
        db.create_unique('server_post_tag_set', ['post_id', 'tag_id'])

        # Adding model 'Blog'
        db.create_table('server_blog', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('author', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('url', self.gf('django.db.models.fields.URLField')(max_length=500)),
        ))
        db.send_create_signal('server', ['Blog'])

        # Adding model 'Related'
        db.create_table('server_related', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('source', self.gf('django.db.models.fields.related.ForeignKey')(related_name='source', to=orm['server.Post'])),
            ('target', self.gf('django.db.models.fields.related.ForeignKey')(related_name='target', to=orm['server.Post'])),
        ))
        db.send_create_signal('server', ['Related'])

        # Adding model 'Visit'
        db.create_table('server_visit', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('ip', self.gf('django.db.models.fields.GenericIPAddressField')(default='', max_length=39, null=True, blank=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal('server', ['Visit'])

        # Adding model 'View'
        db.create_table('server_view', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(related_name='user_views', to=orm['auth.User'])),
            ('post', self.gf('django.db.models.fields.related.ForeignKey')(related_name='post_views', to=orm['server.Post'])),
            ('date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal('server', ['View'])

        # Adding model 'PostBody'
        db.create_table('server_postbody', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('post', self.gf('django.db.models.fields.related.ForeignKey')(related_name='bodies', to=orm['server.Post'])),
            ('content', self.gf('django.db.models.fields.TextField')(max_length=10000)),
            ('html', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal('server', ['PostBody'])

        # Adding model 'PostRevision'
        db.create_table('server_postrevision', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('post', self.gf('django.db.models.fields.related.ForeignKey')(related_name='revisions', to=orm['server.Post'])),
            ('diff', self.gf('django.db.models.fields.TextField')()),
            ('content', self.gf('django.db.models.fields.TextField')()),
            ('author', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal('server', ['PostRevision'])

        # Adding model 'Note'
        db.create_table('server_note', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sender', self.gf('django.db.models.fields.related.ForeignKey')(related_name='note_sender', to=orm['auth.User'])),
            ('target', self.gf('django.db.models.fields.related.ForeignKey')(related_name='note_target', to=orm['auth.User'])),
            ('content', self.gf('django.db.models.fields.CharField')(default='', max_length=5000)),
            ('html', self.gf('django.db.models.fields.CharField')(default='', max_length=5000)),
            ('date', self.gf('django.db.models.fields.DateTimeField')(db_index=True)),
            ('unread', self.gf('django.db.models.fields.BooleanField')(default=True, db_index=True)),
            ('type', self.gf('django.db.models.fields.IntegerField')(default=1)),
            ('url', self.gf('django.db.models.fields.URLField')(default='', max_length=200, blank=True)),
        ))
        db.send_create_signal('server', ['Note'])

        # Adding model 'Vote'
        db.create_table('server_vote', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('author', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('post', self.gf('django.db.models.fields.related.ForeignKey')(related_name='votes', to=orm['server.Post'])),
            ('type', self.gf('django.db.models.fields.IntegerField')(db_index=True)),
            ('date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, db_index=True, blank=True)),
        ))
        db.send_create_signal('server', ['Vote'])

        # Adding model 'Badge'
        db.create_table('server_badge', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('type', self.gf('django.db.models.fields.IntegerField')()),
            ('unique', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('secret', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('count', self.gf('django.db.models.fields.IntegerField')(default=0)),
        ))
        db.send_create_signal('server', ['Badge'])

        # Adding model 'Award'
        db.create_table('server_award', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('badge', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['server.Badge'])),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('date', self.gf('django.db.models.fields.DateTimeField')()),
        ))
        db.send_create_signal('server', ['Award'])


    def backwards(self, orm):
        # Deleting model 'UserProfile'
        db.delete_table('server_userprofile')

        # Deleting model 'Tag'
        db.delete_table('server_tag')

        # Deleting model 'Post'
        db.delete_table('server_post')

        # Removing M2M table for field tag_set on 'Post'
        db.delete_table('server_post_tag_set')

        # Deleting model 'Blog'
        db.delete_table('server_blog')

        # Deleting model 'Related'
        db.delete_table('server_related')

        # Deleting model 'Visit'
        db.delete_table('server_visit')

        # Deleting model 'View'
        db.delete_table('server_view')

        # Deleting model 'PostBody'
        db.delete_table('server_postbody')

        # Deleting model 'PostRevision'
        db.delete_table('server_postrevision')

        # Deleting model 'Note'
        db.delete_table('server_note')

        # Deleting model 'Vote'
        db.delete_table('server_vote')

        # Deleting model 'Badge'
        db.delete_table('server_badge')

        # Deleting model 'Award'
        db.delete_table('server_award')


    models = {
        'auth.group': {
            'Meta': {'object_name': 'Group'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        'auth.permission': {
            'Meta': {'ordering': "('content_type__app_label', 'content_type__model', 'codename')", 'unique_together': "(('content_type', 'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['contenttypes.ContentType']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.Group']", 'symmetrical': 'False', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        'server.award': {
            'Meta': {'object_name': 'Award'},
            'badge': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['server.Badge']"}),
            'date': ('django.db.models.fields.DateTimeField', [], {}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['auth.User']"})
        },
        'server.badge': {
            'Meta': {'object_name': 'Badge'},
            'count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'secret': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'type': ('django.db.models.fields.IntegerField', [], {}),
            'unique': ('django.db.models.fields.BooleanField', [], {'default': 'False'})
        },
        'server.blog': {
            'Meta': {'object_name': 'Blog'},
            'author': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['auth.User']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'url': ('django.db.models.fields.URLField', [], {'max_length': '500'})
        },
        'server.note': {
            'Meta': {'object_name': 'Note'},
            'content': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '5000'}),
            'date': ('django.db.models.fields.DateTimeField', [], {'db_index': 'True'}),
            'html': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '5000'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sender': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'note_sender'", 'to': "orm['auth.User']"}),
            'target': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'note_target'", 'to': "orm['auth.User']"}),
            'type': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'unread': ('django.db.models.fields.BooleanField', [], {'default': 'True', 'db_index': 'True'}),
            'url': ('django.db.models.fields.URLField', [], {'default': "''", 'max_length': '200', 'blank': 'True'})
        },
        'server.post': {
            'Meta': {'object_name': 'Post'},
            'accepted': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'answer_count': ('django.db.models.fields.IntegerField', [], {'default': '0', 'blank': 'True'}),
            'author': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['auth.User']"}),
            'changed': ('django.db.models.fields.BooleanField', [], {'default': 'True', 'db_index': 'True'}),
            'content': ('django.db.models.fields.TextField', [], {'max_length': '10000'}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'db_index': 'True'}),
            'full_score': ('django.db.models.fields.IntegerField', [], {'default': '0', 'db_index': 'True', 'blank': 'True'}),
            'html': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'lastedit_date': ('django.db.models.fields.DateTimeField', [], {}),
            'lastedit_user': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'editor'", 'to': "orm['auth.User']"}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'children'", 'null': 'True', 'to': "orm['server.Post']"}),
            'rank': ('django.db.models.fields.FloatField', [], {'default': '0', 'blank': 'True'}),
            'root': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'descendants'", 'null': 'True', 'to': "orm['server.Post']"}),
            'score': ('django.db.models.fields.IntegerField', [], {'default': '0', 'db_index': 'True', 'blank': 'True'}),
            'slug': ('django.db.models.fields.SlugField', [], {'max_length': '200', 'blank': 'True'}),
            'status': ('django.db.models.fields.IntegerField', [], {'default': '100'}),
            'tag_set': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['server.Tag']", 'symmetrical': 'False', 'blank': 'True'}),
            'tag_val': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'title': ('django.db.models.fields.TextField', [], {'max_length': '200'}),
            'type': ('django.db.models.fields.IntegerField', [], {'db_index': 'True'}),
            'url': ('django.db.models.fields.URLField', [], {'default': "''", 'max_length': '200', 'blank': 'True'}),
            'views': ('django.db.models.fields.IntegerField', [], {'default': '0', 'db_index': 'True', 'blank': 'True'})
        },
        'server.postbody': {
            'Meta': {'object_name': 'PostBody'},
            'content': ('django.db.models.fields.TextField', [], {'max_length': '10000'}),
            'html': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'post': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'bodies'", 'to': "orm['server.Post']"})
        },
        'server.postrevision': {
            'Meta': {'object_name': 'PostRevision'},
            'author': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['auth.User']"}),
            'content': ('django.db.models.fields.TextField', [], {}),
            'date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'diff': ('django.db.models.fields.TextField', [], {}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'post': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'revisions'", 'to': "orm['server.Post']"})
        },
        'server.related': {
            'Meta': {'object_name': 'Related'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'source': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'source'", 'to': "orm['server.Post']"}),
            'target': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'target'", 'to': "orm['server.Post']"})
        },
        'server.tag': {
            'Meta': {'object_name': 'Tag'},
            'count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {'max_length': '50', 'db_index': 'True'})
        },
        'server.userprofile': {
            'Meta': {'object_name': 'UserProfile'},
            'about_me': ('django.db.models.fields.TextField', [], {'default': "''", 'null': 'True', 'blank': 'True'}),
            'about_me_html': ('django.db.models.fields.TextField', [], {'default': "''", 'null': 'True', 'blank': 'True'}),
            'bronze_badges': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'display_name': ('django.db.models.fields.CharField', [], {'default': "'User'", 'max_length': '250', 'db_index': 'True'}),
            'gold_badges': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'last_visited': ('django.db.models.fields.DateTimeField', [], {}),
            'location': ('django.db.models.fields.TextField', [], {'default': "''", 'null': 'True', 'blank': 'True'}),
            'my_tags': ('django.db.models.fields.TextField', [], {'default': "''", 'max_length': '250', 'null': 'True', 'blank': 'True'}),
            'new_messages': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'scholar': ('django.db.models.fields.TextField', [], {'default': "''", 'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'score': ('django.db.models.fields.IntegerField', [], {'default': '0', 'blank': 'True'}),
            'silver_badges': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'status': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'type': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "'profile'", 'unique': 'True', 'to': "orm['auth.User']"}),
            'uuid': ('django.db.models.fields.TextField', [], {'unique': 'True', 'db_index': 'True'}),
            'website': ('django.db.models.fields.URLField', [], {'default': "''", 'max_length': '250', 'null': 'True', 'blank': 'True'})
        },
        'server.view': {
            'Meta': {'object_name': 'View'},
            'date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'post': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'post_views'", 'to': "orm['server.Post']"}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'user_views'", 'to': "orm['auth.User']"})
        },
        'server.visit': {
            'Meta': {'object_name': 'Visit'},
            'date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'ip': ('django.db.models.fields.GenericIPAddressField', [], {'default': "''", 'max_length': '39', 'null': 'True', 'blank': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['auth.User']"})
        },
        'server.vote': {
            'Meta': {'object_name': 'Vote'},
            'author': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['auth.User']"}),
            'date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'db_index': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'post': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'votes'", 'to': "orm['server.Post']"}),
            'type': ('django.db.models.fields.IntegerField', [], {'db_index': 'True'})
        }
    }

    complete_apps = ['server']