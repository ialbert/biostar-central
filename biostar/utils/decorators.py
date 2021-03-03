import logging, functools, time
from functools import partial
from django.conf import settings
from django.shortcuts import redirect
from django.contrib import messages
import sys

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


def dtimer():
    """
    Return a disabled timer.
    """

    class inner(object):
        def __init__(self, secs, **kwargs):
            self.secs = secs

        def __call__(self, f, *args, **kwargs):
            pass

    return inner


def btimer():
    """
    Return blocking timer, runs the function once.
    """

    class inner(object):
        def __init__(self, secs, **kwargs):
            self.secs = secs

        def __call__(self, func, *args, **kwargs):
            func(*args, **kwargs)

    return inner


def ttimer():
    """
    Return threaded timer
    """

    class inner(object):
        def __init__(self, secs, **kwargs):
            self.secs = secs

        def __call__(self, func, *args, **kwargs):
            # The loop repeats the timer.
            def loop():
                ticker = threading.Event()
                while not ticker.wait(self.secs):
                    func(*args, **kwargs)

            # Run process in separate thread, once.
            logger.info(f"new time thread for function f{func} {args} {kwargs}")
            t = threading.Thread(target=loop, daemon=True)
            t.start()

    return inner


def utimer():
    """
    Ensure uwsgi is installed and return the timer
    """
    from uwsgidecorators import timer

    return timer


def ctimer():
    """
    Ensure celery is installed and construct a callable object to match interface.
    Inside, we dynamically add this function to the beat schedule.

    Adopted from:
    https://docs.celeryproject.org/en/master/userguide/periodic-tasks.html#beat-entries
    """

    from biostar.celery import app

    print("CALLED HERE ONCE")

    class inner(object):

        def __init__(self, secs, **kwargs):
            self.secs = secs

        # Handler means that not to evaluate the app at module level when calling f()
        # on_after_finalize
        @app.on_after_configure.connect
        def __call__(self, f, *args, **kwargs):
            f = app.task(f)
            # Add the entry to the beat_schedule setting.
            app.add_periodic_task(schedule=self.secs, sig=f, kwargs=kwargs, args=args,
                                  name=f.__name__)

    return inner


def thread(*args, **kwargs):
    def outer(func, **kwargs):
        @functools.wraps(func)
        def inner(*args, **kwargs):
            # Run process in separate thread.
            logger.info(f"new thread for function f{func} {args} {kwargs}")
            t = threading.Thread(target=func, args=args, kwargs=kwargs, daemon=True)
            t.start()

        # Gains an attribute called spool that runs the function in the background.
        inner.spool = inner
        inner.delay = inner
        return inner

    return outer


def spooler():
    """
    Return a uwsgi spooler after ensuring
    """
    # Ensure uwsgi is installed.

    from uwsgidecorators import spool

    def inner(f):
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

    return inner


def cworker():
    # Ensure celery is installed
    from biostar.celery import app

    def inner(f):
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
        worker = app.task(f)
        # Compatible with uwsgi interface.
        worker.spool = worker.delay
        return worker

    return inner


def block():
    def outer(func, *args, **kwargs):
        @functools.wraps(func)
        def inner(*args, **kwargs):
            return func(*args, **kwargs)

        inner.spool = inner
        inner.delay = inner

        return inner

    return outer


def disabled():
    def outer(func, *args, **kwargs):
        @functools.wraps(func)
        def inner(*args, **kwargs):
            pass

        inner.spool = inner
        inner.delay = inner
        return inner

    return outer


def threaded():
    def inner(f):
        worker = thread()(f)
        return worker

    return inner


def mappings():
    mapper = {
        'block': {'worker': block, 'timer': btimer},
        'uwsgi': {'worker': spooler, 'timer': utimer},
        'celery': {'worker': cworker, 'timer': ctimer},
        'threaded': {'worker': threaded, 'timer': ttimer},
        'disabled': {'worker': disabled, 'timer': dtimer},
    }

    if settings.TASK_RUNNER not in mapper:
        logger.error(f"Invalid Task. valid options : {mapper.keys()}")
        raise Exception('Invalid task.')

    logger.info(f'tasks and timers set to {settings.TASK_RUNNER}')

    return mapper


# def pick_task():
#     """
#     Return decorator used to run tasks, uses settings.TASK_RUNNER.
#     """
#     # Validate and retrieve map of workers.
#     opts = mappings()
#     # Call primary function once here and return worker decorator.
#     picked = opts.get(settings.TASK_RUNNER)[0]()
#     return picked
#

def retrieve(ftype):
    # Validate and retrieve map of workers.
    opts = mappings()
    picked = opts.get(settings.TASK_RUNNER)[ftype]()
    return picked

#
# def pick_timer():
#     """
#
#     Return decorator used to run timer, uses settings.TASK_RUNNER.
#     """
#     # Validate and retrieve map of workers.
#     opts = mappings()
#     # Call primary function once here and return worker decorator.
#     picked = opts.get(settings.TASK_RUNNER)[1]()
#     return picked


try:
    # Pick what task and timer functions to use.
    SPOOLER = retrieve('worker')
    TIMER = retrieve('timer')

except Exception as exc:
    # Disable tasks when there are errors
    SPOOLER = disabled()
    TIMER = dtimer()
    logger.error(f'Error picking task: {settings.TASK_RUNNER}, {exc}. Tasks disabled.')


def task(f):
    """
    Select task decorator depending on settings.TASK_RUNNER.
    """
    return SPOOLER(f)


def timer(f):
    """
    Select timer decorator depending on settings.TASK_RUNNER.
    """
    return TIMER(f)
