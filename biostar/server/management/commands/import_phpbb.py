# import_phpbb v1.0
# This script imports the a PHPBB based forum into the Biostars instance
#
# 
# The script uses the backup of the database in 'csv' format
# Required tables - 'phpbb_posts' & 'phpbb_users'
# 
#
# Example Backup - MySQL:
#	mysql -B -u $USERNAME -p"$PASSWORD" $DBNAME  --execute="SELECT * FROM $TABLENAMEW;" |sed "s/'/\'/;s/\t/\",\"/g;s/^/\"/;s/$/\"/;s/\n//g" > $TABLENAME.csv
#
#
# Script Usage:
#
#	python manage.py import_phpbb -p /path/to/posts/table.csv -u /path/to/users/table.csv
#
#	Example: pyhton manage.py import_phpbb -p ./phpbb_posts.csv -u ./phpbb_users.csv
#

from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
import sys, logging, os, datetime, pytz
from django.core.exceptions import ImproperlyConfigured
from datetime import date
from django.utils import timezone
from django.conf import settings
import csv, re
from biostar.apps.users.models import User
from biostar.apps.posts.models import Post

logger = logging.getLogger('simple-logger')

class Command(BaseCommand):
    help = 'Import phpbb posts from csv backup\nUsage: python manage.py import_phpbb -p /path/to/posts/table.csv -u /path/to/users/table.csv'

    option_list = BaseCommand.option_list + (
        make_option("-p", '--post', dest='post', default=False, help='path to posts table csv'),
        make_option("-u", '--user', dest='user', default=False, help='path to users table csv'),
    )

    def handle(self, *args, **options):
        fname = options['post']
        uname = options['user']

        if fname and uname:
            if fname.endswith('.csv') and uname.endswith('.csv'):
                import_posts(fname,uname)
            else:
                logger.error('Wrong format! Please provide a csv file (.csv)')
        else:
            if not fname:
                logger.error('No posts file name supplied')
            if not uname:
                logger.error('No users file name supplied')
            logger.info('try -h for more help')


#Returns the filteres body content
def bfilter(body):
    body = body.replace('\\n','\n')
    body = re.sub(r'\[url:\w+\]','URL: ',body)
    body = re.sub(r'\[url=','URL: ',body)
    body = re.sub(r'\[/url:\w+\]','',body)
    return body


#Extracts and returns all the posts
def get_posts(fname):
    f = open(fname,'rb')
    rows = csv.reader(f)
    rows.next()

    allposts = []
    for row in rows:
        uid = row[3]
        title = row[14]
        body = bfilter(row[15])
        date = datetime.datetime.fromtimestamp(int(row[6]), tz=pytz.utc)
        temp = [uid,title,body,date]
        allposts.append(temp)

    f.close()
    return allposts


#Returns all users from file
def get_all_users(uname):
    f = open(uname,'rb')
    rows = csv.reader(f)
    rows.next()

    allusers = []
    for row in rows:
        uid = int(row[0])
        name = row[7]
        email = row[12]
        temp = [uid,name,email]
        allusers.append(temp)

    f.close()
    return allusers


#Extracts user information with user_id
def get_user(uid, allusers):

    found=0
    for row in allusers:
        if row[0] == uid:
            name = row[1]
            email = row[2]
            found=1
            break

    try:
        user = User.objects.get(email=email)
        logger.info('Fetched user: %s' % user.name)
    except:
        user = User(email=email, name=name)
        user.save()
        logger.info('Created user: %s' % user.name)

    return user 


#Imports the posts into django models
def import_posts(fname, uname):

    logger.info('Extracting posts from file...')
    allposts = get_posts(fname)

    logger.info('Extracting users from file...')
    allusers = get_all_users(uname)

    post_count=0
    for single in allposts:
        uid = int(single[0])
        title = single[1]
        body = single[2]
        date = single[3]
        logger.info('Fetched post : %s' % title)

        #Getting the post user
        user = get_user(uid,allusers)

        #List to hold answer posts which could not be matched
        orphans = []

        if title.startswith('Re:'):
            try:
                ptitle = title[4:]
                parent = Post.objects.get(title=ptitle)
                post = Post(title=title, content=body, author=user, type= Post.ANSWER, creation_date=date)
                post.parent=parent
                post.root=parent
                post.save()
                post_count+=1
            except:
                post = Post(title=title, author=user, type= Post.ANSWER, content=body, creation_date=date)
                orphans.append(post)
        else:
            post = Post(title=title, content=body, author=user, type = Post.QUESTION, creation_date=date)
            post.save()
            post_count+=1


    #Now try to match posts which could not be matched before
    if orphans:
        logger.info('Matching posts which could not be matched earlier')
        for post in orphans:
            try:
                title = post.title
                ptitle = title[4:]
                parent = Post.objects.get(title__startswith=ptitle)
                post.parent=parent
                post.root=parent
                post.save()
                post_count+=1
            except:
                pass

    print post_count, ' posts created'
    logger.info('DONE!')


if __name__ == '__main__':
    pass