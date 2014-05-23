# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Peer'
        db.create_table(u'tracker_peer', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('info_hash', self.gf('django.db.models.fields.CharField')(max_length=40, db_index=True)),
            ('peer_id', self.gf('django.db.models.fields.CharField')(max_length=40, db_index=True)),
            ('user_agent', self.gf('django.db.models.fields.CharField')(max_length=120, blank=True)),
            ('ip_address', self.gf('django.db.models.fields.IPAddressField')(max_length=15)),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=40, blank=True)),
            ('port', self.gf('django.db.models.fields.IntegerField')()),
            ('uploaded', self.gf('django.db.models.fields.PositiveIntegerField')(default=0)),
            ('downloaded', self.gf('django.db.models.fields.PositiveIntegerField')(default=0)),
            ('left', self.gf('django.db.models.fields.PositiveIntegerField')(default=0)),
            ('lastupdate_date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
            ('creation_date', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['users.User'], null=True)),
        ))
        db.send_create_signal(u'tracker', ['Peer'])


    def backwards(self, orm):
        # Deleting model 'Peer'
        db.delete_table(u'tracker_peer')


    models = {
        u'sites.site': {
            'Meta': {'ordering': "(u'domain',)", 'object_name': 'Site', 'db_table': "u'django_site'"},
            'domain': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'tracker.peer': {
            'Meta': {'object_name': 'Peer'},
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'downloaded': ('django.db.models.fields.PositiveIntegerField', [], {'default': '0'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'info_hash': ('django.db.models.fields.CharField', [], {'max_length': '40', 'db_index': 'True'}),
            'ip_address': ('django.db.models.fields.IPAddressField', [], {'max_length': '15'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '40', 'blank': 'True'}),
            'lastupdate_date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'left': ('django.db.models.fields.PositiveIntegerField', [], {'default': '0'}),
            'peer_id': ('django.db.models.fields.CharField', [], {'max_length': '40', 'db_index': 'True'}),
            'port': ('django.db.models.fields.IntegerField', [], {}),
            'uploaded': ('django.db.models.fields.PositiveIntegerField', [], {'default': '0'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['users.User']", 'null': 'True'}),
            'user_agent': ('django.db.models.fields.CharField', [], {'max_length': '120', 'blank': 'True'})
        },
        u'users.user': {
            'Meta': {'object_name': 'User'},
            'activity': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'badges': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'email': ('django.db.models.fields.EmailField', [], {'unique': 'True', 'max_length': '255', 'db_index': 'True'}),
            'flair': ('django.db.models.fields.CharField', [], {'default': "u''", 'max_length': '15'}),
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

    complete_apps = ['tracker']