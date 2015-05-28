from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
import sys, logging, os
from django.core.exceptions import ImproperlyConfigured
from datetime import date
from django.utils import timezone
from django.conf import settings
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


#Fetches the environment variable
def get_env(name, func=None):
    """Get the environment variable or return exception"""
    try:
        if func:
            return func(os.environ[name])
        else:
            return unicode(os.environ[name], encoding="utf-8")
    except KeyError:
        msg = "*** Required environment variable %s not set." % name
        raise ImproperlyConfigured(msg)


#Imports the posts into django models
def import_posts(allposts):
	from biostar.apps.users.models import User
	from biostar.apps.posts.models import Post

	emailhost=get_env('EMAIL_HOST')
	email = 'sqlimport@' + emailhost + '.com'
	try:
		u = User.objects.get(email=email)
	except:
		u = User(email=email, name='sqlimport')
		u.save()
		u.profile.date_joined = timezone.now()
		u.profile.last_login = timezone.now()
		u.profile.save()

	post_count=0
	for single in allposts:
		title = single[0]
		body = single[1]
		logger.info('Fetched post : %s' % title)
		if title.startswith('Re:'):
			ptitle = title[4:]
			try:
				parent = Post.objects.get(title=ptitle)
				post = Post(title=title, content=body, author=u, type= Post.ANSWER)
				post.parent=parent
				post.root=parent
				post.save()
				post_count+=1
			except:
				pass
		else:
			post = Post(title=title, content=body, author=u, type = Post.QUESTION)
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