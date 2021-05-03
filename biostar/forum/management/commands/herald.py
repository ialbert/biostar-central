import logging
from typing import Any
import os, sys
from django.template import loader
from django.core.management.base import BaseCommand
from biostar.forum.models import Post, Herald
from biostar.forum.util import now
from biostar.accounts.models import User
from biostar.planet.models import BlogPost, Blog
from biostar.forum import auth
from mistune import Markdown
from django.shortcuts import reverse
from django.conf import settings

logger = logging.getLogger('engine')


def create_blog(heralds):
    """
    Link a group of herald publications to a blog post.
    """
    date = now()
    # Create the Blog post title and content first displayed.
    template = "herald/publication.md"
    tmpl = loader.get_template(template_name=template)
    context = dict(heralds=heralds)
    content = tmpl.render(context)

    # Get the Herald blog where all publications belong to.
    hblog = Blog.objects.filter(link=reverse('herald_list')).first()
    if not hblog:
        logger.warning(f"Herald blog does not exist.")
        return

    # Convert template to html
    html = Markdown()(text=content)

    # Create a blog post
    title = f"Biostar Herald: {date.date()} issue."
    blgpost = BlogPost.objects.create(title=title, blog=hblog, content=content, insert_date=date, html=html,
                                      creation_date=date)
    # Update this blog post link to be the herald issue page.
    blgpost.link = reverse('herald_issue', kwargs=dict(blog_pk=blgpost.pk))
    blgpost.save()

    # Link the herald publications to this blog post.
    heralds.update(blog_post=blgpost)
    user = User.objects.filter(is_superuser=True).first()
    auth.db_logger(user=user, text=f"published {heralds.count()} submissions in {title}")

    return


def run_publisher(limit=20):
    """
    Publish most recently accepted herald_list submissions.
    """

    # Get most recent heralds
    heralds = Herald.objects.filter(status=Herald.ACCEPTED)

    create_blog(heralds)

    # Update the heralds
    heralds.update(status=Herald.PUBLISHED)
    return


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--publish', action='store_true', default=False,
                            help="Create a publication out of the most recently accepted herald_list submissions")
        parser.add_argument('--limit', type=int,
                            help="How many submission to collate in a publication")

    def handle(self, *args, **options):
        # Index all un-indexed posts that have a root.
        logger.debug(f"Database: {settings.DATABASE_NAME}")

        publish = options['publish']

        if publish:
            run_publisher()
            return
