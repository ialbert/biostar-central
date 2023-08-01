"""
Saves a thread to a file.
"""
import json, time
import logging
import os
from itertools import *
from django.conf import settings
from django.contrib import sitemaps
from django.contrib.sitemaps import GenericSitemap
from django.contrib.sites.models import Site
from django.core.management.base import BaseCommand, CommandError
from django.template import loader
from django.utils.encoding import smart_str
from django.core import serializers
import tempfile
from biostar.forum.models import Post, JsonPost

import json
from tqdm import tqdm

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

def print_json(data):
    text = json.dumps(data, indent=4)

    print(text)

def serialize_post(uid):

    # Annotate post with all the ids that have the same root id.
    post = Post.objects.select_related("author__profile", "root", "lastedit_user__profile").get(uid=uid)

    # Select all posts with root id = post.root_id but exlude the post itself.
    root_ids = Post.objects.filter(root=post).exclude(pk=post.pk).values_list('id', flat=True)
    root_ids = list(root_ids)

    #print (desc)

    tag_val = post.tag_val.replace(",", " ")
    tag_val = tag_val.split()

    #fname = create_path(post)
    text = serializers.serialize('json', [ post, ])

    data = json.loads(text)

    # Add a field to the data
    data[0]['fields']['root_ids'] = root_ids
    data[0]['fields']['tag_val'] = tag_val

    return data[0]

def batched(it, n ):
    batch = list(islice(it, n))
    while batch:
        yield batch
        batch = list(islice(it, n))

def save_posts(uid):

    N = 100
    jp = JsonPost.objects.all().first()

    print (jp)



    # Clear the database.
    JsonPost.objects.all().delete()

    total = Post.objects.count()

    posts = Post.objects.all()

    pbar = tqdm(total=total)

    def as_json(post):
        data = serialize_post(uid=post.uid)
        jp = JsonPost(pk=post.pk, data=data)
        return jp

    posts = map(as_json, posts)

    chunks = batched(posts, N)
    count = 0

    for index, chunk in enumerate(chunks):
        JsonPost.objects.bulk_create(chunk, ignore_conflicts=True)
        count += len(chunk)
        pbar.update(len(chunk))
        pbar.set_description("# chunks %s" % count)
        time.sleep(1)

    return




class Command(BaseCommand):
    help = 'Creates a sitemap in the export folder of the site'

    def add_arguments(self, parser):
        parser.add_argument('--uid', default=9164, help="post id")

    def handle(self, *args, **options):
        pid = int(options['uid'])
        save_posts(uid=pid)
