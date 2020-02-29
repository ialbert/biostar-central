import logging, functools

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
                # Run process in separate thread.
                t = threading.Thread(target=func, args=args, kwargs=kwargs, daemon=True)
                t.start()
            inner.spool = inner
            return inner
        # Gains an attribute called spool that
        # falls back to the same function
        return outer

    # Create a synchronous version of the timer
    def timer(secs, **kwargs):
        def outer(func):
            @functools.wraps(func)
            def inner(*args, **kwargs):
                result = func(*args, **kwargs)
                return result
            inner.timer = inner
            return inner
        # Gains an attribute called timer that
        # falls back to the same function
        return outer
