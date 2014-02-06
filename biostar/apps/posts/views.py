# Create your views here.
from django.shortcuts import render_to_response
from django.views.generic import TemplateView, DetailView, ListView, FormView, UpdateView
from .models import Post
from django import forms
from django.core.urlresolvers import reverse
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Fieldset, Submit, ButtonHolder
from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.contrib import messages
from . import auth
from braces.views import LoginRequiredMixin

# Create your views here.
class PostEditForm(forms.Form):
    title = forms.CharField()
    content = forms.CharField(widget=forms.Textarea, required=False)

    def __init__(self, *args, **kwargs):
        super(PostEditForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.layout = Layout(
            Fieldset(
                'Post',
                'title',
                'content',
            ),
            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )



class EditPost(LoginRequiredMixin, FormView):
    """
    Edits a post.
    """

    # The template_name attribute must be specified in the calling apps.
    template_name = ""
    form_class = PostEditForm
    user_fields = "name email".split()
    prof_fields = "location website info scholar".split()

    def get(self, request, *args, **kwargs):
        initial = {}
        pk = int(self.kwargs['pk'])
        if pk > 0:
            post = Post.objects.get(pk=pk)
            initial = dict(title=post.title, content=post.content)
            print (post.id, post.root)

        form = self.form_class(initial=initial)
        return render(request, self.template_name, {'form': form})

    def post(self, request, *args, **kwargs):

        pk = int(self.kwargs['pk'])
        form = self.form_class(request.POST)

        # Validate the form.
        if not form.is_valid():
            return render(request, self.template_name, {'form': form})

        user = request.user
        data = form.cleaned_data
        # Valid forms start here.
        if pk == 0:
            post = Post(
                title=data['title'], content=data['content'], author=user, type=Post.FORUM,
            )
            post.save()

        else:
            post = Post.objects.get(pk=self.kwargs['pk'])
            post = auth.post_permissions(request=request, post=post)
            if not post.is_editable:
                messages.error(request, "This user may not modify the post")
                return HttpResponseRedirect(reverse("home"))
            post.title = data['title']
            post.content  = data['content']
            post.lastedit_user = user
            post.save()
            messages.success(request, "Post updated")

        return HttpResponseRedirect(post.get_absolute_url())

    def get_success_url(self):
        return reverse("user-details", kwargs=dict(pk=self.kwargs['pk']))

