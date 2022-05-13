"""
Saves a thread to a file.
"""
import json
import logging
import os
from django.conf import settings
from django.contrib import sitemaps
from django.contrib.sitemaps import GenericSitemap
from django.contrib.sites.models import Site
from django.core.management.base import BaseCommand, CommandError
from django.template import loader
from django.utils.encoding import smart_str
from django.core import serializers
import tempfile
from biostar.forum.models import Post

logger = logging.getLogger("engine")

def atomic_save(fname, data):

    # Create the directory if it does not exist.
    path = os.path.dirname(fname)

    if not os.path.exists(path):
        logger.info(f"Creating directory {path}")
        os.makedirs(path, exist_ok=True)

    # Atomic write
    temp_file = tempfile.NamedTemporaryFile(delete=False, dir=path, suffix=".tmp")

    with open(temp_file.name, mode="w") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())

    # Move the file to the final location
    os.replace(temp_file.name, fname)

def create_path(post):
    creation_date = post.creation_date
    year = creation_date.strftime("%Y")
    month = creation_date.strftime("%m")
    day = creation_date.strftime("%d")

    print(year, month, day)

    path = f"posts/{year}/{month}/{day}/{post.root.uid}/{post.uid}.json"

    full_path = os.path.join(settings.JSON_DIR, path)

    return full_path

def save_thread(uid):
    post = Post.objects.select_related("author__profile", "root", "lastedit_user__profile").get(uid=uid)

    fname = create_path(post)

    text = serializers.serialize('json', [ post, ])
    data = json.loads(text)
    text = json.dumps(data, indent=4)

    #print(text)

    logger.info(f"Saving {post.uid} to {fname}")

    atomic_save(fname, data=text)


class Command(BaseCommand):
    help = 'Creates a sitemap in the export folder of the site'

    def add_arguments(self, parser):
        parser.add_argument('--uid', default=0, help="post id")

    def handle(self, *args, **options):
        pid = int(options['uid'])
        save_thread(uid=pid)
