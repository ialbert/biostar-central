"""
Saves a thread to a file.
"""
import os
from django.conf import settings
from django.contrib.sitemaps import GenericSitemap
from django.contrib.sites.models import Site
from biostar.forum.models import Post
from django.utils.encoding import smart_str
from django.template import loader
from django.core.management.base import BaseCommand, CommandError
import logging
from django.contrib import sitemaps
import json

logger = logging.getLogger("engine")

def save_threead(uid):
    post = Post.objects.select_related("author__profile", "lastedit_user__profile").get(uid=uid)

    creation_date = post.creation_date
    year = creation_date.strftime("%Y")
    month = creation_date.strftime("%m")
    day = creation_date.strftime("%d")

    print(year, month, day)

    data = dict(
        id=post.id,
        root_id=post.root_id,
        parent_id=post.parent_id,
        uid=post.uid,
        path = f"{year}/{month}/{day}",
        title=post.title,
        type=post.type,
        type_name = post.get_type_display(),
        creation_date=post.creation_date.isoformat(),
        lastedit_date=post.lastedit_date.isoformat(),
        lastedit_user=post.lastedit_user_id,
        author_id=post.author_id,
        author_name=post.author.profile.name,
        content=post.content,
        html=post.html,
    )

    text = json.dumps(data, indent=4)
    print(text)


class Command(BaseCommand):
    help = 'Creates a sitemap in the export folder of the site'

    def add_arguments(self, parser):
        parser.add_argument('--uid', default=0, help="post id")

    def handle(self, *args, **options):
        pid = int(options['uid'])
        save_threead(uid=pid)