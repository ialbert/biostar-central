# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting field 'Badge.secret'
        db.delete_column(u'badges_badge', 'secret')

        # Deleting field 'Badge.description'
        db.delete_column(u'badges_badge', 'description')

        # Adding field 'Badge.desc'
        db.add_column(u'badges_badge', 'desc',
                      self.gf('django.db.models.fields.CharField')(default='', max_length=200),
                      keep_default=False)

        # Adding field 'Badge.icon'
        db.add_column(u'badges_badge', 'icon',
                      self.gf('django.db.models.fields.CharField')(default='fa fa-asterisk', max_length=250),
                      keep_default=False)


    def backwards(self, orm):
        # Adding field 'Badge.secret'
        db.add_column(u'badges_badge', 'secret',
                      self.gf('django.db.models.fields.BooleanField')(default=False),
                      keep_default=False)

        # Adding field 'Badge.description'
        db.add_column(u'badges_badge', 'description',
                      self.gf('django.db.models.fields.CharField')(default='', max_length=200),
                      keep_default=False)

        # Deleting field 'Badge.desc'
        db.delete_column(u'badges_badge', 'desc')

        # Deleting field 'Badge.icon'
        db.delete_column(u'badges_badge', 'icon')


    models = {
        u'badges.award': {
            'Meta': {'object_name': 'Award'},
            'badge': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['badges.Badge']"}),
            'date': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['users.User']"})
        },
        u'badges.badge': {
            'Meta': {'object_name': 'Badge'},
            'count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'desc': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '200'}),
            'icon': ('django.db.models.fields.CharField', [], {'default': "'fa fa-asterisk'", 'max_length': '250'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'type': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'unique': ('django.db.models.fields.BooleanField', [], {'default': 'False'})
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

    complete_apps = ['badges']