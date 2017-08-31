from django.db import models
from django.contrib.auth.models import User


class Post(models.Model):

    author = models.ForeignKey(User)
    title = models.CharField(max_length=140)
    body = models.TextField()
    date = models.DateTimeField()

    def __str__(self):
        return self.title
    
    
    

