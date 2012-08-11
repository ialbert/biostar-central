"""
Extracts parameters for the biostar.sh script
"""
import sys
from django.conf import settings

if __name__ == '__main__':
    key = sys.argv[1]
    
    engine = settings.DATABASES['default']['ENGINE']
    is_sqlite = engine.endswith('sqlite3')
    is_postgres = engine.endswith('psycopg2')
    
    if key == 'PG_NAME':
        PG_NAME = settings.DATABASES['default']['NAME'] if is_postgres else ""
        print PG_NAME
    elif key == "SQLITE_DBNAME":
        SQLITE_DBNAME = settings.DATABASES['default']['NAME'] if is_sqlite else ""
        print SQLITE_DBNAME
    else:
        print settings.DATABASES['default'][key]