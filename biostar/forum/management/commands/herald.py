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
from django.shortcuts import reverse
from django.conf import settings

logger = logging.getLogger('engine')


def create_blogs(heralds):
    """
    Link a group of herald publications to a blog post.
    """
    date = now()
    # Create the Blog post title and content first displayed.
    title = f"Biostar Herard: {date.date()} issue."
    template = "herald/publication.md"
    tmpl = loader.get_template(template_name=template)
    context = dict(title=title, heralds=heralds)
    content = tmpl.render(context)

    # Get the Herald blog where all publications belong to.
    hblog = Blog.objects.filter(link=reverse('herald_list')).first()
    if not hblog:
        logger.warning(f"Herald blog does not exist.")
        return

    # Create a blog post
    blgpost = BlogPost.objects.create(title=title, blog=hblog, content=content, insert_date=date, creation_date=date)
    # Update this blog post link to be the herald issue page.
    blgpost.link = reverse('herald_issue', kwargs=dict(blog_pk=blgpost.pk))
    blgpost.save()

    # Link the publications to this blog post.
    heralds.update(blog_post=blgpost)
    user = User.objects.filter(is_superuser=True).first()
    auth.db_logger(user=user, text=f"created heard:{title}")

    return


def publish(limit=20):
    """
    Publish most recently accepted herald_list submissions.
    """

    # Get most recent heralds
    heralds = Herald.objects.filter(status=Herald.ACCEPTED)

    # Update the heralds
    heralds.update(status=Herald.PUBLISHED)

    create_blogs(heralds)
    return


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--issue', action='store_true', default=False,
                            help="Create a publication out of the most recently accepted herald_list submissions")
        parser.add_argument('--limit', type=int,
                            help="How many submission to collate in a publication")

    def handle(self, *args, **options):
        # Index all un-indexed posts that have a root.
        logger.debug(f"Database: {settings.DATABASE_NAME}")

        issue = options['issue']

        if issue:
            publish()
            return
