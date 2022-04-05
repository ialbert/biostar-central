import logging, functools, time, os
from functools import wraps
from functools import partial
from ratelimit.decorators import ratelimit
from django.conf import settings
from django.http import Http404
from django.shortcuts import redirect
from django.contrib import messages
import sys

logger = logging.getLogger('engine')
import threading


def is_moderator(f):
    def inner(request, **kwargs):
        user = request.user
        if user.is_authenticated and user.profile.is_moderator:
            return f(request, **kwargs)
        messages.warning(request, "You need to be a moderator to perform this action.")
        return redirect('/')

    return inner



def authenticated(func):

    def _wrapper_(request, **kwargs):
        if request.user.is_anonymous:
            messages.error(request, "You need to be logged in to view this page.")
            return redirect("/")

        return func(request, **kwargs)

    return _wrapper_

def is_staff(f):
    """
    Only run functions with the
    """
    def inner(request, **kwargs):
        user = request.user
        if user.is_authenticated and (user.is_staff or user.is_superuser):
            return f(request, **kwargs)
        messages.warning(request, "You need to be a staff member to perform this action.")
        return redirect('/')

    return inner


def reset_count(key):
    """
    Set value of given key in settings.SESSION_COUNT_KEY to 0.
    """
    session_key = settings.SESSION_COUNT_KEY

    def outer(func):
        @wraps(func)
        def inner(request, **kwargs):

            if request.user.is_authenticated:
                # Get the count from sessions
                counts = request.session.get(session_key, {})
                counts[key] = 0
                request.session[session_key] = counts

            results = func(request, **kwargs)

            return results

        return inner

    return outer


def check_params(allowed):
    """
    Validate if only allowed params are present in request.GET.
    """
    def outer(func):
        @wraps(func)
        def inner(request, **kwargs):
            incoming = set(request.GET.keys())
            # Expected parameter names.
            diff = incoming - allowed
            if diff:
                logger.warning(f"invalid get request parameters {diff}")
                raise Http404("Parameter does not exist.")

            return func(request, **kwargs)

        return inner

    return outer


def limited(key, rate):
    """
    Make a blocking rate limiter that does not raise an exception
    """
    def outer(func):

        @ratelimit(key=key, rate=rate)
        def inner(request, **kwargs):

            was_limited = getattr(request, 'limited', False)
            if was_limited:
                msg = "Too many requests from same IP address. Temporary ban."
                messages.warning(request, msg)
                raise Http404(msg)

            return func(request, **kwargs)

        return inner

    return outer


def timeit(func):
    """
    Print how long function takes.
    """

    @wraps(func)
    def inner(*args, **kwargs):
        start = time.time()
        val = func(*args, **kwargs)

        delta = int((time.time() - start) * 1000)
        msg = f"time={delta}ms for {func.__name__}"

        if delta > 1000:
            msg = f'SLOW: {msg}'
            logger.info(msg)
        else:
            logger.debug(msg)

        return val

    return inner


def check_lock(lock):
    """
    Check if lock directory exists before calling function
    """

    def __inner(func):

        def __wrapper(*args, **kwargs):

            if os.path.isdir(lock):
                logger.warning('Lock directory detected, function is already running')
                sys.exit()

            # Try to run function
            try:
                # Make the lock directory
                os.makedirs(lock, exist_ok=True)
                out = func(*args, **kwargs)
            except Exception as exc:
                logger.error(exc)
                out = None
            finally:
                # Clean the locks.
                os.rmdir(lock)

            # Return function output
            return out

        return __wrapper

    return __inner


def d_timer():
    """
    Return disabled timer.
    """

    class inner(object):
        def __init__(self, secs, **kwargs):
            self.secs = secs

        def __call__(self, f, *args, **kwargs):
            pass

    return inner


def b_timer():
    """
    Return blocking timer
    """

    class inner(object):
        def __init__(self, secs, **kwargs):
            self.secs = secs

        def __call__(self, func, *args, **kwargs):
            func(*args, **kwargs)

    return inner


def t_timer():
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


def u_timer():
    """
    Return uwsgi timer
    """
    from uwsgidecorators import timer

    return timer


def c_timer():
    """
    Construct a celery timer decorator.
    Inside the __call__, it dynamically adds the given function to the beat schedule.

    Adopted from:
    https://docs.celeryproject.org/en/master/userguide/periodic-tasks.html#beat-entries
    """

    from biostar.celery import app

    class inner(object):

        def __init__(self, secs, **kwargs):
            self.secs = secs

        # Handler means that not to evaluate the app at module level when calling f()
        # on_after_finalize
        @app.on_after_configure.connect
        def __call__(self, f, *args, **kwargs):
            f = app.task(f)
            # Add the entry to the beat_schedule.
            app.add_periodic_task(schedule=self.secs,
                                  sig=f,
                                  kwargs=kwargs,
                                  args=args,
                                  name=f.__name__)

    return inner


def thread(*args, **kwargs):
    """
    Return a threaded worker
    """
    def outer(func, **kwargs):
        @functools.wraps(func)
        def inner(*args, **kwargs):
            # Run process in separate thread.
            logger.debug(f"new thread for function f{func} {args} {kwargs}")
            t = threading.Thread(target=func, args=args, kwargs=kwargs, daemon=True)
            t.start()

        # Gains an attribute called spool that runs the function in the background.
        inner.spool = inner
        inner.delay = inner
        return inner

    return outer


def u_worker():
    """
    Return a uwsgi spooler compatible with celery interface
    """
    # Ensure uwsgi is installed.

    from uwsgidecorators import spool

    def inner(f):
        """
        Alias to call .spool when calling .delay

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


def c_worker():
    """
    Return a celery worker compatible with uwsgi interface
    """
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


def b_worker():
    """
    Return a blocking decorator that runs the function once.
    """
    def outer(func, *args, **kwargs):
        @functools.wraps(func)
        def inner(*args, **kwargs):
            logger.debug(f"running f{func} {args} {kwargs}")
            return func(*args, **kwargs)

        inner.spool = inner
        inner.delay = inner

        return inner

    return outer


def d_worker():
    """
    Return a d_worker decorator that does nothing
    """
    def outer(func, *args, **kwargs):
        @functools.wraps(func)
        def inner(*args, **kwargs):
            pass

        inner.spool = inner
        inner.delay = inner
        return inner

    return outer


def t_worker():
    """
    Wrap a threaded worker and to match interface.
    """
    def inner(f):
        worker = thread()(f)
        return worker

    return inner


def select_runner(name):
    """
    Return runner based on name ( worker or timer ) and settings.TASK_RUNNER.
    """
    mapper = {
        'block': {'worker': b_worker, 'timer': b_timer},
        'uwsgi': {'worker': u_worker, 'timer': u_timer},
        'celery': {'worker': c_worker, 'timer': c_timer},
        'threaded': {'worker': t_worker, 'timer': t_timer},
        'disable': {'worker': d_worker, 'timer': d_timer},
    }

    if settings.TASK_RUNNER not in mapper:
        logger.error(f"Invalid Task. valid options : {mapper.keys()}")
        raise Exception('Invalid task.')

    # Call primary function here and return worker decorator.
    decorator = mapper.get(settings.TASK_RUNNER)[name]()
    return decorator


try:
    # Initiate the runners
    WORKER = select_runner('worker')
    TIMER = select_runner('timer')
    logger.debug(f'workers and timers set to {settings.TASK_RUNNER}')

except Exception as exc:
    # Disable tasks when there are errors, raising exceptions breaks migration.
    WORKER = d_worker()
    TIMER = d_timer()
    logger.warning(f'Error initializing task: {settings.TASK_RUNNER}.')
    logger.warning(f'Tasks disabled: {exc}.')


def task(f):
    """
    Utility function to access worker decorator.
    """
    return WORKER(f)


def timer(f):
    """
    Utility function to access timer decorator.
    """
    return TIMER(f)
