from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
import sys, logging, os
from django.core.exceptions import ImproperlyConfigured
from datetime import date
from django.utils import timezone
from django.conf import settings
from biostar3.forum.models import *
import csv


logger = logging.getLogger('simple-logger')

class Command(BaseCommand):
    help = 'Import from .sql backup'

    option_list = BaseCommand.option_list + (
        make_option("-f", '--file', dest='file', default=False, help='import file'),
    )

    def handle(self, *args, **options):
        fname = options['file']

        if fname:
        	if fname.endswith('.csv'):
        		import_bb_file(fname)
        	else:
        		logger.error('Wrong format! Please provide a csv file (.csv)')
        else:
        	if not fname:
        		logger.error('No file name supplied')
        	logger.info('try -h for more help')


#Extracts and returns posts from the csv
def extract(fname):
	f = open(fname, 'rb')
	rows = csv.reader(f)

	allposts = []
	for row in rows:
		title = row[14]
		body = row[15]
		temp=[title,body]
		allposts.append(temp)

	allposts.pop(0)
	return allposts


#Imports the posts into django models
def import_posts(allposts):

    # Get the default group.
    default_group = UserGroup.objects.get(domain=settings.DEFAULT_GROUP_DOMAIN)

    # Get the first admin user
    admin = User.objects.get(email=settings.ADMINS[0][1])

    post_count=0
    for single in allposts:
        title = single[0]
        body = single[1]
        logger.info('Fetched post : %s' % title)
        if title.startswith('Re:'):
            ptitle = title[4:]
            try:
                parent = Post.objects.get(title=ptitle)
                post = Post(title=title, author=admin, type= Post.ANSWER, content=body, usergroup=default_group)
                post.parent=parent
                post.root=parent
                post.save()
                post_count+=1
            except:
                pass
        else:
            post = Post(title=title, author=admin, type = Post.QUESTION, content=body, usergroup=default_group)
            post.save()
            post_count+=1
    logger.info('%d posts created' % post_count)    


#Governs the improt process
def import_bb_file(fname):
	logger.info('Extracting posts from file...')
	allposts = extract(fname)

	logger.info('Importing posts...')
	import_posts(allposts)
	logger.info('DONE!')


if __name__ == '__main__':
    pass