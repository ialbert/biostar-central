import re
import bleach
import logging
import time
import uuid
from itertools import islice, count
from datetime import datetime
from django.utils.timezone import utc


def fixcase(name):
    return name.upper() if len(name) == 1 else name.lower()


def now():
    return datetime.utcnow().replace(tzinfo=utc)


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def strip_tags(text):
    "Strip html tags from text"
    text = bleach.clean(text, tags=[], attributes={}, styles=[], strip=True)
    return text


def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)


def timer_func():
    """
    Prints progress on inserting elements.
    """

    last = time.time()

    def elapsed(msg):
        nonlocal last
        now = time.time()
        sec = round(now - last, 1)
        last = now
        print(f"{msg} in {sec} seconds")

    def progress(index, step=5000, msg=""):
        nonlocal last
        if index % step == 0:
            print("**" * 5)
            print()
            elapsed(f"... {index} {msg}")
            print()
            print("**" * 5)

    return elapsed, progress