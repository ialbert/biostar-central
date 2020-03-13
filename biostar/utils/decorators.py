import logging, functools
from django.conf import settings
logger = logging.getLogger('biostar')
import threading

try:
    # Loads up uwsgi
    from uwsgidecorators import spool, timer

except Exception as exc:
    logger.warning("uwsgi decorators not found, tasks are synchronous")

    # Create a synchronous version of the spooler
    def spool(pass_arguments=True):
        def outer(func):
            @functools.wraps(func)
            def inner(*args, **kwargs):

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

    # Create a synchronous version of the timer
    def timer(secs, **kwargs):
        def outer(func):
            @functools.wraps(func)
            def inner(*args, **kwargs):
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
