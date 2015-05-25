from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
import sys, logging, os
import MySQLdb as mdb 
from django.core.exceptions import ImproperlyConfigured
from datetime import date
from django.utils import timezone
from django.conf import settings

logger = logging.getLogger('simple-logger')

class Command(BaseCommand):
    help = 'Import from .sql backup'

    option_list = BaseCommand.option_list + (
        make_option("-f", '--file', dest='file', default=False, help='import file'),
        make_option("-l", '--host', dest='host', default='localhost', help='mysql hostname'),
        make_option("-u", '--user', dest='user', default='root', help='mysql user'),
        make_option("-p", '--pass', dest='pass', default=None, help='mysql password (will not be saved)'),
    )

    def handle(self, *args, **options):
        fname = options['file']
        host = options['host']
        user = options['user']
        password = options['pass']

        if fname and password:
            import_sql_file(fname, host, user, password)
        else:
        	if not fname:
        		logger.error('No file name supplied')
        	if not password:
        		logger.error('MySQL password needed')
        	logger.info('try -h for more help')


#Removes unwanted lines from sql file
def sanitize(filename):
	f = open(filename,'r')
	sqltemp=[]
	for i in f.readlines():
	    if i[0]!='-' and i[0]!='/' and i[0]!='\n':
	    	if i[-1:] == '\n':
	    		i=i[:-1]
	        sqltemp.append(i)
	sql=[]
	temp=[]
	for i in sqltemp:
		if i[-1:]==';':
			temp.append(i)
			temp = ''.join(temp)
			sql.append(temp)
			temp=[]
		else:
			temp.append(i)
	return sql


#Creates/Delets the temporary mysql database
def temp_db(host, user, password, cmd):
	if cmd == 'c':
		sql = 'CREATE DATABASE import_temp'
	if cmd == 'd':
		sql = 'DROP DATABASE import_temp'
	conn = mdb.connect(host,user,password)
	cursor = conn.cursor()
	cursor.execute(sql)
	conn.commit()
	conn.close()


#Imports the sql file to the temporary database
def import_to_temp(sql, host, user, password):
	conn = mdb.connect(host,user,password,'import_temp')
	cursor = conn.cursor()
	cursor._defer_warnings = True
	for i in sql:
		cursor.execute(i)
	conn.commit()
	conn.close()


#Imports the posts into django models
def import_posts(host, user, password):
	from biostar3.forum.models import *


	# Get the default group.
	default_group = UserGroup.objects.get(domain=settings.DEFAULT_GROUP_DOMAIN)

	# Get the first admin user
	admin = User.objects.get(email=settings.ADMINS[0][1])

	conn = mdb.connect(host,user,password,'import_temp')
	cursor = conn.cursor()
	cursor._defer_warnings = True
	sql = 'select * from phpbb_posts;'
	try:
		cursor.execute(sql)
		results = cursor.fetchall()
	except:
		logger.error('Unable to fetch posts from temp_db')
		sys.exit()

	post_count=0
	for result in results:
		title = result[14]
		body = result[15]
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

	conn.close()


#Governs the import process
def import_sql_file(filename, host, user, password):
	temp_db(host, user, password, 'c')
	logger.info('Created temporary database')

	logger.info('Sanitizing input file')
	sql = sanitize(filename)

	logger.info('Importing to temp db...')
	import_to_temp(sql, host, user, password)

	logger.info('Importing posts...')
	import_posts(host, user, password)

	temp_db(host, user, password, 'd')
	logger.info('Deleted temporary database')


if __name__ == '__main__':
    pass