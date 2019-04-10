import logging, time, shutil, subprocess
from django.core import management
from django.utils.encoding import force_text


import time

logger = logging.getLogger("engine")

HAS_UWSGI = False


COUNTER = 1

try:
    from uwsgidecorators import *

    HAS_UWSGI = True


    def post_created(user):
        "Executed on a post creation"
        logger.info("post created")


    @app.task
    def call_command(name, *args, **kwargs):
        "Calls a django command in a delayed fashion"
        logger.info("calling django command %s with %s and %s" % (name, args, kwargs))
        from django.core.management import call_command
        call_command(name, *args, **kwargs)


    @app.task
    def test(*args, **kwds):
        logger.info("*** executing task %s %s, %s" % (__name__, args, kwds))



except ModuleNotFoundError as exc:
    pass