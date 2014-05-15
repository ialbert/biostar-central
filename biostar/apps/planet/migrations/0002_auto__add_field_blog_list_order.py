# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'Blog.list_order'
        db.add_column(u'planet_blog', 'list_order',
                      self.gf('django.db.models.fields.IntegerField')(default=0),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting field 'Blog.list_order'
        db.delete_column(u'planet_blog', 'list_order')


    models = {
        u'planet.blog': {
            'Meta': {'object_name': 'Blog'},
            'active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'desc': ('django.db.models.fields.TextField', [], {'default': "''", 'blank': 'True'}),
            'feed': ('django.db.models.fields.URLField', [], {'max_length': '200'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'link': ('django.db.models.fields.URLField', [], {'max_length': '200'}),
            'list_order': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'title': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '255'})
        },
        u'planet.blogpost': {
            'Meta': {'object_name': 'BlogPost'},
            'blog': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['planet.Blog']"}),
            'content': ('django.db.models.fields.TextField', [], {'default': "''", 'max_length': '20000'}),
            'creation_date': ('django.db.models.fields.DateTimeField', [], {'db_index': 'True'}),
            'html': ('django.db.models.fields.TextField', [], {'default': "''"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'insert_date': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'db_index': 'True'}),
            'link': ('django.db.models.fields.URLField', [], {'max_length': '200'}),
            'published': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'uid': ('django.db.models.fields.CharField', [], {'max_length': '200'})
        }
    }

    complete_apps = ['planet']