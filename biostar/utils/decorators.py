import logging, functools
from django.conf import settings
logger = logging.getLogger('biostar')
import threading

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
    def spool(pass_arguments=True):
        def outer(func):
            @functools.wraps(func)
            def inner(*args, **kwargs):
                if settings.DISABLE_TASKS:
                    return
                if settings.MULTI_THREAD:
                    # Run process in separate thread.
                    logger.info(f"new thread for function f{func} {args} {kwargs}")
                    t = threading.Thread(target=func, args=args, kwargs=kwargs, daemon=True)
                    t.start()
                else:
                    func(*args, **kwargs)
            inner.spool = inner
            return inner
        # Gains an attribute called spool that runs the function in the background.
        return outer

    # Create a threaded version of the timer
    def timer(secs, **kwargs):
        def outer(func):
            @functools.wraps(func)
            def inner(*args, **kwargs):
                if settings.DISABLE_TASKS:
                    return

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

            inner.timer = inner
            return inner
        # Gains an attribute called timer that will run the function periodically.
        return outer
