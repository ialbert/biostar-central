"""
Extracts parameters for the biostar.sh script
"""
import sys
from django.conf import settings


if __name__ == '__main__':
    
    if len(sys.argv)!= 2:
        print ('(!) one value may be extracted at a time' )
        sys.exit()

    key = sys.argv[1]
    engine = settings.DATABASES['default']['ENGINE']
    is_postgres = engine.endswith('psycopg2')

    if is_postgres:
        remap = dict(PG_DBNAME="NAME", PG_USER="USER", PG_PASSWD="PASSWORD")
    else:
        remap = dict(SQLITE_DBNAME="NAME")
    
    key = remap.get(key, '')
    value =  settings.DATABASES['default'].get(key, '')
    
    print value