try:
    from .celery import app as celery_app
    __all__ = ['celery_app']
except ImportError as exc:
    print(f'celery error {exc}')

VERSION = '2.3.5'
