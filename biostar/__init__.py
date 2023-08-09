try:
    from .celery import app as celery_app
    __all__ = ['celery_app']
except ImportError as exc:
    pass

VERSION = '2.3.6'

# Monkeypatch the translation module to avoid errors in Django 4.2
import django.utils.translation
django.utils.translation.ugettext_lazy = django.utils.translation.gettext_lazy

import django.utils.encoding

#django.utils.encoding.force_text = django.utils.encoding.force_str
#django.utils.encoding.smart_text = django.utils.encoding.smart_str
