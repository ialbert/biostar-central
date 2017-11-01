from biostar.engine.settings import *

DEBUG=True

WSGI_APPLICATION = 'conf.devel.devel_wsgi.application'

ALLOWED_HOSTS = [ 'localhost', 'www.lvh.me' ]
