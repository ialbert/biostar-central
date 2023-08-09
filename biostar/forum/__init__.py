# Monkeypatch the translation module to avoid errors in Django 4.2
import django.utils.translation
django.utils.translation.ugettext_lazy = django.utils.translation.gettext_lazy

