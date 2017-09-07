from django.db import models
from django.contrib.auth.models import User
    
#class Title(models.Model):
#    title = models.CharField(max_length=140)
#
#    def __str__(self):
#        return self.title


class Project(models.Model):

    title = models.CharField(max_length=256)
    author = models
    text = models.TextField()
    html = models.TextField()
    date = models.DateTimeField()

    def __str__(self):
        return self.title


class Data(models.Model):
    title = models.CharField(max_length=256)
    author = models
    text = models.TextField()
    html = models.TextField()
    date = models.DateTimeField()

    def __str__(self):
        return self.title

