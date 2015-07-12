"""
Crawls over a directory and all subdirectories and adds pages that it finds
to the current django instance.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import os, logging, sys, re
from .compat import *
from django.db import transaction
from biostar3.forum import auth
from biostar3.forum.models import FlatPage, User, Post
from django.template.loader import get_template
from django.template import Context, Template

from django.conf import settings

logger = logging.getLogger('biostar')

# Match Django template comments.
PATTERN = re.compile(r'^{#\s+(?P<name>\w+)\s?=\s?(?P<value>[\S\s]+) #}')


def parse_metadata(raw_content):
    "Attempts to parse out metadata from django comments"
    lines = raw_content.splitlines()[:20]
    meta = dict()
    lines = map(strip, lines)
    for line in lines:

        m = PATTERN.search(line)
        if m:
            name, value = m.group('name'), m.group('value')
            meta[name] = value
    return meta

def render_page(raw_content, params={}):
    return Template(raw_content).render(Context(params))

def add_one(user, raw_content, path, update=False):
    errors = []

    if not raw_content:
        errors.append("no content in {}".format(path))

    meta = parse_metadata(raw_content)
    slug = meta.get("slug")
    if not slug:
        errors.append("slug field is missing from {}".format(path))
    title = meta.get("title")
    if not title:
        errors.append("title field is missing from {}".format(path))

    # Check for the slug.
    page = FlatPage.objects.filter(slug=slug).first()
    if page and (not update):
        errors.append("slug {} already exists from {}".format(slug, path))

    if errors:
        return errors
    else:
        with transaction.atomic():
            logger.info("creating: {}".format(slug))
            rendered_content = render_page(raw_content)
            data = dict(
                title=title, type=Post.PAGE, content=rendered_content
            )
            post = auth.create_toplevel_post(data=data, user=user)
            page = FlatPage.objects.create(post=post, slug=slug)
            return []

def add_all(path, update=False):
    valid_exts = {".html", ".md"}

    admin = User.objects.filter(email=settings.ADMINS[0][1]).first()

    for dirpath, dirnames, filenames in os.walk(path):
        for name in sorted(filenames):
            start, ext = os.path.splitext(name)
            if ext not in valid_exts:
                continue

            path = os.path.join(dirpath, name)
            content = open(path).read()

            errors = add_one(admin, content, path)

            for error in errors:
                logger.error(error)

if __name__ == '__main__':
    add_all(sys.argv[1])
