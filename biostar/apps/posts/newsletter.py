__author__ = 'ialbert'

from django.views.generic import TemplateView, DetailView, ListView, FormView, UpdateView, CreateView
from .models import Post, EmailSub, EmailEntry
from braces.views import LoginRequiredMixin, StaffuserRequiredMixin
from django.http import HttpResponseRedirect, Http404


class EntryList(LoginRequiredMixin, StaffuserRequiredMixin, ListView):
    model = EmailEntry
    template_name = "newsletter/entry_list.html"

    def get_context_data(self, **kwargs):
        context = super(EntryList, self).get_context_data(**kwargs)
        return context


class CreateEntry(LoginRequiredMixin, StaffuserRequiredMixin, CreateView):
    model = EmailEntry
    template_name = "newsletter/create_entry.html"


class UdpdateEntry(LoginRequiredMixin, UpdateView):
    model = EmailEntry
    template_name = "newsletter/create_entry.html"