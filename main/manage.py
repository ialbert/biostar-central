#!/usr/bin/env python
from django.core.management import execute_manager
import imp
try:
    imp.find_module('settings') # Assumed to be in the same directory.
except ImportError:
    import sys
    sys.stderr.write("Error: Can't find the file 'settings.py' in the directory containing %r. It appears you've customized things.\nYou'll have to run django-admin.py, passing it your settings module.\n" % __file__)
    sys.exit(1)

import settings

# monkey patching Django to allow make the debug server multithreaded
# as seen in https://code.djangoproject.com/ticket/3357
# new versions of Django will be multithreaded by default
import SocketServer
import django.core.servers.basehttp
"""
django.core.servers.basehttp.WSGIServer = \
    type('WSGIServer',
         (SocketServer.ThreadingMixIn,
          django.core.servers.basehttp.WSGIServer,
          object),
         {})
"""

if __name__ == "__main__":
    execute_manager(settings)
