import logging
from django import forms
from django.template import loader
from django.conf import settings
from biostar.accounts.models import User, Profile
from django.contrib import messages
from django.shortcuts import render, redirect, reverse
from biostar.planet.models import Blog, BlogPost
from django.db.models import F
from biostar.forum import auth, util
from biostar.forum.models import Post, SharedLink
from biostar.emailer.models import EmailGroup, EmailSubscription
from biostar.utils.decorators import is_moderator, authenticated
from biostar.forum.tasks import herald_emails

from .const import *

logger = logging.getLogger("engine")

MIN_CHARS = 5
MAX_CONTENT = 15000
MIN_CONTENT = 5
MAX_TITLE = 400
MAX_TAGS = 5


class HeraldSubmit(forms.Form):
    url = forms.URLField(required=True)
    text = forms.CharField(widget=forms.Textarea(attrs=dict(rows='5')), max_length=MAX_CONTENT, required=False,
                           strip=False)

    # Maximum amount of submissions.
    MAX = 5

    def __init__(self, user=None, *args, **kwargs):
        self.user = user
        super(HeraldSubmit, self).__init__(*args, **kwargs)

    def clean(self):
        cleaned_data = super(HeraldSubmit, self).clean()

        if self.user.is_anonymous:
            raise forms.ValidationError("You need to be logged in.")

        # if exists:
        #    raise forms.ValidationError("This link already exists.")

        # Check if non mode user is over max submissions
        count = SharedLink.objects.filter(author=self.user, status=SharedLink.SUBMITTED).count()

        if not self.user.profile.is_moderator and count >= self.MAX:
            raise forms.ValidationError(
                f"You already have {count} links submitted, please wait until some are accepted.")

        return cleaned_data


def render_template(template, context):
    tmpl = loader.get_template(template_name=template)
    content = tmpl.render(context)
    return content


def herald_blog(post):
    """
    Create a herald blog post from a post.

    """

    # Get the Biostar Herald blog.
    hlink = reverse('post_topic', kwargs=dict(topic='herald'))
    blog = Blog.objects.filter(link=hlink).first()

    content = post.content.split('.')[:2]
    content = '.'.join(content) + '...'

    # Create the blog post
    BlogPost.objects.create(title=post.title, blog=blog, link=post.get_absolute_url(),
                            uid=post.uid, content=content,
                            creation_date=post.creation_date, insert_date='', published=True)
    return


def remove_declined(limit=50):
    """
    Remove rejected post herald
    """

    declined = SharedLink.objects.filter(status=SharedLink.DECLINED)[:limit]

    # Can not apply .delete() to sliced query set
    SharedLink.objects.filter(pk__in=declined.values_list('pk')).delete()


def herald_publisher(request, limit=20, nmin=1):
    """
    Create one publication from Herald accepted submissions ( up to 'limit' ).
    """

    # Reset status on published links.
    # SharedLink.objects.filter(status=SharedLink.PUBLISHED).update(status=SharedLink.ACCEPTED)

    heralds = SharedLink.objects.filter(status=SharedLink.ACCEPTED).order_by('-pk')[:limit]
    count = heralds.count()

    if count < nmin:
        logger.warning(f"Not enough stories to publish, minimum of {nmin} required.")
        return

    # Create herald content
    date = util.now()

    date_fmt = date.strftime("%A, %B %d, %Y")

    title = f"The Biostar Herald for {date_fmt}"

    port = f':{settings.HTTP_PORT}' if settings.HTTP_PORT else ''

    base_url = f"{settings.PROTOCOL}://{settings.SITE_DOMAIN}{port}"
    authors = set(h.author for h in heralds)
    editors = set(h.editor for h in heralds)

    subscribe_url = reverse('herald_subscribe')
    context = dict(heralds=heralds, title=title, site_domain=settings.SITE_DOMAIN, protocol=settings.PROTOCOL,
                   base_url=base_url, authors=authors, editors=editors, subscribe_url=subscribe_url)

    content = render_template(template="herald/herald_content.md", context=context)

    # Create herald post
    user = User.objects.filter(is_superuser=True).first()
    post = auth.create_post(title=title, content=content, author=user, tag_val='herald', ptype=Post.HERALD,
                            nodups=False, request=request)

    # Tie these submissions to herald post
    hpks = heralds.values_list('pk', flat=True)
    SharedLink.objects.filter(pk__in=hpks).update(status=SharedLink.PUBLISHED, post=post, lastedit_date=date)

    # Log the action
    auth.db_logger(user=user, text=f"published {count} submissions in {title}")

    # Create a herald blog post
    herald_blog(post=post)

    # Send out herald emails
    herald_emails.spool(uid=post.uid)

    # Clean up declined links
    remove_declined()

    return post


@authenticated
def herald_list(request):
    """
    List latest herald_list items
    """

    # List newly submitted links.
    show = request.GET.get('show', 'upcoming')
    # Filter for or exclude published links
    if show == 'published':
        stories = SharedLink.objects.filter(status=SharedLink.PUBLISHED)
    else:
        stories = SharedLink.objects.exclude(status=SharedLink.PUBLISHED)

    stories = stories.order_by('-creation_date')
    stories = stories.select_related('author', 'author__profile')
    user = request.user
    form = HeraldSubmit(user=user)

    if request.method == 'POST':

        form = HeraldSubmit(data=request.POST, user=user)

        if form.is_valid():
            # Add the Link attribute.
            link = form.cleaned_data['url']
            text = form.cleaned_data['text']

            # Create the herald_list objects.
            herald = SharedLink.objects.create(author=user, text=text, url=link)
            messages.success(request, 'Submitted for review.')
            return redirect(reverse('herald_list'))

    # Apply a limit
    stories = stories[:settings.HERALD_LIST_COUNT]
    count = SharedLink.objects.filter(author=user, status=SharedLink.SUBMITTED).count()
    allow_submit = count < form.MAX
    context = dict(stories=stories, tab='herald_list', show=show, form=form, allow_submit=allow_submit)

    return render(request, 'herald/herald_base.html', context)


@is_moderator
def herald_publish(request):
    post = herald_publisher(request)

    if not post:
        messages.error(request, "Not enough submissions to publish.")
        return redirect(reverse('herald_list'))

    return redirect(reverse('post_view', kwargs=dict(uid=post.uid)))
