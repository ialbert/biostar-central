import logging, functools
from functools import partial
from django.conf import settings
from django.shortcuts import redirect
from django.contrib import messages
logger = logging.getLogger('biostar')
import threading


def is_moderator(f):

    def inner(request, **kwargs):
        user = request.user
        if user.is_authenticated and user.profile.is_moderator:
            return f(request, **kwargs)
        messages.warning(request, "You need to be a moderator to perform this action.")
        return redirect('/')

    return inner


def thread(*args, **kwargs):
    def outer(func, **kwargs):
        if settings.DISABLE_TASKS:
            return

        @functools.wraps(func)
        def inner(*args, **kwargs):
            if settings.MULTI_THREAD:
                # Run process in separate thread.
                logger.info(f"new thread for function f{func} {args} {kwargs}")
                t = threading.Thread(target=func, args=args, kwargs=kwargs, daemon=True)
                t.start()
            else:
                func(*args, **kwargs)

        @functools.wraps(func)
        def timed(secs, **kwargs):
            # The loop repeats the timer.
            def loop():
                ticker = threading.Event()
                while not ticker.wait(secs):
                    func(*args, **kwargs)

            if settings.MULTI_THREAD:
                # Run process in separate thread, once.
                logger.info(f"new time thread for function f{func} {args} {kwargs}")
                t = threading.Thread(target=loop, daemon=True)
                t.start()
            else:
                func(*args, **kwargs)

        # Gains an attribute called spool that runs the function in the background.
        inner.spool = inner
        inner.delay = inner
        inner.timer = timed
        return inner

    return outer


try:
    # When run with uwsgi the tasks will be spooled via uwsgi.
    from uwsgidecorators import spool, timer

except Exception as exc:
    #
    # With no uwsgi module the tasks will be spooled.
    # Creating threaded versions of the decorators from uwsgi.
    #
    logger.warning("uwsgi module not found, tasks will run in threads")

    # Create a threaded version of the spooler
    spool = thread
    timer = thread


def spooler(f):
    """
    Alias to call .delay when calling .spool

    @spooler
    def job(foo):
       pass

    # Uwsgi type of launch
    job.spool(foo='')

    # Celery type of launch
    job.delay(foo='')

    """
    worker = spool(pass_arguments=True)(f)
    # Compatible with celery interface.
    worker.delay = worker.spool
    return worker


def threaded(f):
    worker = thread()(f)
    return worker

