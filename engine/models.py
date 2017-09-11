from django.db import models
from django.contrib.auth.models import User
from django.utils import timezone

def make_html(text):
    return text

class Base(models.Model):

    title = models.CharField(max_length=256)
    owner = models.ForeignKey(User)
    text = models.TextField(default='text')
    html = models.TextField(default='html')
    date = models.DateTimeField()

    def __str__(self):
        return self.title

    def save(self, *args, **kwargs):
        now = timezone.now()
        self.date = self.date or now
        self.html = make_html(self.text)
        super(Base, self).save(*args, **kwargs)

    class Meta:
        abstract = True

        
class Data(Base):

    def save(self, *args, **kwargs):
        super(Data, self).save(*args, **kwargs)


class Analysis(Base):

    def save(self, *args, **kwargs):
        super(Analysis, self).save(*args, **kwargs)


class Project(Base):

    data = models.ManyToManyField(Data)
    analysis = models.ManyToManyField(Analysis)

    def save(self, *args, **kwargs):
        super(Project, self).save(*args, **kwargs)






    
