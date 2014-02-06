# Create your views here.
from django.shortcuts import render_to_response
from django.views.generic import TemplateView, DetailView, ListView, FormView, UpdateView
from .models import User
from . import auth
from django import forms
from django.core.urlresolvers import reverse
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit

class UserEditForm(forms.Form):
    name = forms.CharField()
    location = forms.CharField()
    website = forms.CharField()
    scholar = forms.CharField()

    info = forms.CharField(widget=forms.Textarea)

    # form.helper.form_action = reverse('user-edit', args=[event.id])
    # form.helper.form_action = reverse('url_name', kwargs={'book_id': book.id})

    def send_email(self):
        # send email using the self.cleaned_data dictionary
        pass

    def __init__(self, *args, **kwargs):
        super(UserEditForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.add_input(Submit('submit', 'Submit'))


class EditUser(FormView):
    """
    Edits a user.
    """

    form_class = UserEditForm
    # the template name must be specified in the calling apps.
    # template_name = "???"


    """
    def get_object(self):
        obj = super(UserDetails, self).get_object()
        obj = auth.user_permissions(user=self.request.user, target=obj)
        return obj

    def get_context_data(self, **kwargs):
        context = super(UserDetails, self).get_context_data(**kwargs)
        return context
    """