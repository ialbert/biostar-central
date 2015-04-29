"""
Crawls over a directory and all subdirectories and adds pages that it finds
to the current django instance.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import os, logging, sys, re
from .compat import *
from django.db import transaction
from biostar3.forum import auth
from biostar3.forum.models import UserGroup, FlatPage, User, Post
from django.template.loader import get_template
from django.template import Context, Template

from django.conf import settings

logger = logging.getLogger('biostar')

# Match Django template comments.
PATTERN = re.compile(r'^{#\s+(?P<name>\w+)\s?=\s?(?P<value>[\S\s]+) #}')


def parse_metadata(path):
    "Attempts to parse out metadata from django comments"
    lines = open(path).read().splitlines()[:20]
    meta = dict()
    lines = map(strip, lines)
    for line in lines:

        m = PATTERN.search(line)
        if m:
            name, value = m.group('name'), m.group('value')
            meta[name] = value
    return meta

def render_page(path, params={}):
    content = open(path).read()
    templ = Template(content)
    context = Context(params)
    content = templ.render(context)
    return content

def add_all(path, update=False):
    valid_exts = {".html", ".md"}

    usergroup = UserGroup.objects.filter(domain=settings.DEFAULT_GROUP_DOMAIN).first()

    admin = User.objects.filter(email=settings.ADMINS[0][1]).first()

    for dirpath, dirnames, filenames in os.walk(path):
        for name in sorted(filenames):
            start, ext = os.path.splitext(name)
            if ext not in valid_exts:
                continue

            path = os.path.join(dirpath, name)

            meta = parse_metadata(path)
            slug = meta.get("slug")
            title = meta.get("title")

            if not slug:
                logger.error("slug field is missing from {}".format(path))
                continue

            if not title:
                logger.error("title field is missing from {}".format(path))
                continue


            # Check for the slug.
            page = FlatPage.objects.filter(slug=slug).first()

            if not page or update:
                with transaction.atomic():
                    logger.info("creating: {}".format(slug))
                    content = render_page(path)
                    data = dict(
                        title=title, type=Post.PAGE, content=content
                    )
                    post = auth.create_toplevel_post(data=data, user=admin, group=usergroup)
                    page = FlatPage.objects.create(post=post, slug=slug)


if __name__ == '__main__':
    add_all(sys.argv[1])
