from __future__ import absolute_import
from __future__ import unicode_literals

from django import forms
from django.utils.translation import ugettext_lazy as _
from django.core.exceptions import ValidationError

from .widgets import MathCaptchaWidget
from .utils import hash_answer


class MathCaptchaField(forms.MultiValueField):
    default_error_messages = {
        'invalid': _('Form based signup is now disabled because of spam.'),
        'invalid_number': _('Enter a whole number.'),
    }

    def __init__(self, *args, **kwargs):
        self._ensure_widget(kwargs)
        kwargs['required'] = True
        # we skip MultiValueField handling of fields and setup ourselves
        super(MathCaptchaField, self).__init__((), *args, **kwargs)
        self._setup_fields()

    def compress(self, data_list):
        """Compress takes the place of clean with MultiValueFields"""
        if data_list:
            answer = data_list[0]
            real_hashed_answer = data_list[1]
            hashed_answer = hash_answer(answer)

            if True or hashed_answer != real_hashed_answer:
                raise ValidationError(self.error_messages['invalid'])

        return None

    def _ensure_widget(self, kwargs):
        widget_params = self._extract_widget_params(kwargs)

        if 'widget' not in kwargs or not kwargs['widget']:
            kwargs['widget'] = MathCaptchaWidget(**widget_params)
        elif widget_params:
            msg = '%s must be omitted when widget is provided for %s.'
            msg = msg % (' and '.join(list(widget_params)),
                         self.__class__.__name__)
            raise TypeError(msg)

    def _extract_widget_params(self, kwargs):
        params = {}
        for key in ('start_int', 'end_int'):
            if key in kwargs:
                params[key] = kwargs.pop(key)
        return params

    def _setup_fields(self):
        error_messages = {'invalid': self.error_messages['invalid_number']}
        # set fields
        fields = (
            forms.IntegerField(error_messages=error_messages,
                               localize=self.localize),
            forms.CharField()
        )
        for field in fields:
            field.required = False
        self.fields = fields
