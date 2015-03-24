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

def get_group(meta):
    domain=meta.get("domain")
    usergroup = UserGroup.objects.filter(domain=domain).first()
    if not usergroup:
        logger.info("default usergroup")
        usergroup = UserGroup.objects.filter(domain=settings.DEFAULT_GROUP_DOMAIN).first()
    return usergroup


def crawl(path, update=False):
    valid_exts = { ".html" }

    admin = User.objects.filter(email=settings.ADMINS[0][1]).first()

    for dirpath, dirnames, filenames in os.walk(path):
        for name in sorted(filenames):
            start, ext = os.path.splitext(name)
            if ext not in valid_exts:

                continue

            path = os.path.join(dirpath, name)
            logger.info("processing {}".format(path))
            meta = parse_metadata(path)
            slug = meta.get("slug")
            title = meta.get("title")
            usergroup = get_group(meta)

            if not slug:
                logger.error("slug field is missing from {}".format(path))
            if not title:
                logger.error("title field is missing from {}".format(path))

            FlatPage.objects.filter(post__usergroup=usergroup, slug=slug).delete()

            page = FlatPage.objects.filter(post__usergroup=usergroup, slug=slug).first()
            if not page:
                with transaction.atomic():
                    logger.info("creating {}".format(slug))
                    content = open(path).read()
                    templ = Template(content)
                    params = dict()
                    context = Context(params)
                    content = templ.render(context)
                    data = dict(
                        title=title, type=Post.PAGE, content=content
                    )

                    post = auth.create_toplevel_post(data=data, user=admin,group=usergroup)
                    page = FlatPage.objects.create(post=post, slug=slug)

                continue




if __name__ == '__main__':
    crawl(sys.argv[1])
