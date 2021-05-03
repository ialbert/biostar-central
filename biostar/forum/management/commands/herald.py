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


def run_publisher(limit=20):
    """
    Create one publication from Herald accepted submissions ( up to 'limit' ).
    """

    # Get most recent heralds
    heralds = Herald.objects.filter(status=Herald.ACCEPTED)[:limit]

    date = now()
    # Create content used when listing blog posts.
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

    # Create an issue of the herald
    title = f"Biostar Herald: {date.date()} "
    blgpost = BlogPost.objects.create(title=title, blog=hblog, content=content, insert_date=date, html=html,
                                      creation_date=date)
    # Blog post link points to a herald issue url.
    blgpost.link = reverse('herald_issue', kwargs=dict(blog_pk=blgpost.pk))
    blgpost.save()

    # Link the herald publications to this blog post.
    heralds.update(blog_post=blgpost, status=Herald.PUBLISHED)
    # Log the action.
    user = User.objects.filter(is_superuser=True).first()
    auth.db_logger(user=user, text=f"published {heralds.count()} submissions in {title}")


class Command(BaseCommand):
    help = 'Create search index for the forum app.'

    def add_arguments(self, parser):
        parser.add_argument('--publish', action='store_true', default=False,
                            help="Create a publication out of the most recently accepted herald_list submissions")
        parser.add_argument('--limit', type=int,
                            help="How many submission to collate in this publication", default=20)

    def handle(self, *args, **options):
        # Index all un-indexed posts that have a root.
        logger.debug(f"Database: {settings.DATABASE_NAME}")

        publish = options['publish']
        limit = options['limit']
        if publish:
            run_publisher(limit=limit)
            return
