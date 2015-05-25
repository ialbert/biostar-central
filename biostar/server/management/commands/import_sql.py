from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
import sys, logging
import MySQLdb as mdb 

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
            init_import(fname, host, user, password)
        else:
        	if not fname:
        		logger.error('No file name supplied')
        	if not password:
        		logger.error('MySQL password needed')
        	logger.info('try -h for more help')


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

def import_to_temp(sql, host, user, password):
	conn = mdb.connect(host,user,password,'import_temp')
	cursor = conn.cursor()
	cursor._defer_warnings = True
	for i in sql:
		cursor.execute(i)
	conn.commit()
	conn.close()

def import_sql_file(filename, host, user, password):
	temp_db(host, user, password, 'c')
	logger.info('Created temporary database')

	logger.info('Sanitizing input file')
	sql = sanitize(filename)

	logger.info('Importing to temp db...')
	import_to_temp(sql, host, user, password)

	temp_db(host, user, password, 'd')
	logger.info('Deleted temporary database')

def init_import(fname, host, user, password):
	import_sql_file(fname, host, user, password)


if __name__ == '__main__':
    pass