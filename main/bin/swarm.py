__author__ = 'ialbert'
from django.conf import settings
from datetime import datetime, timedelta
from main.server import const
from functools import *
import operator
import time
from main.server import models

SKIP = "skip"

def time_stamp(date):
    tp = date.timetuple()
    ts = int(time.mktime(tp)) * 1000
    return ts

def encode(text):
    text = text.encode("ascii", errors="replace")
    return text

class Event(object):
    "An object that encapsulates the event related information"
    xml_head  = '<?xml version="1.0"?>'
    start_tag = '<file_events>'
    end_tag   = '</file_events>'
    xml_patt  = '<event filename="%(label)s" date="%(time_stamp)s" author="%(display_name)s" type="%(type)s" on="%(nice_date)s"/>'

    def __init__(self, post):
        # core attributes
        self.allow = True
        self.selected = False
        self.post = post
        self.author = post.author
        self.date = post.creation_date
        self.type = 0

        # derived attributes
        self.nice_date = self.date.strftime("%b %d, %Y")
        self.label = "/%s/%s" % (self.post.get_type_display(), self.post.id)
        self.time_stamp = time_stamp(self.date)
        self.display_name = encode(self.author.profile.display_name)

    def xml(self):
        return Event.xml_patt % self.__dict__

def fetch(limit=10):
    "generates the XML"
    posts = models.Post.objects.all().select_related("author__profile").exclude(type=const.POST_BLOG).order_by("creation_date")
    posts = posts[:limit]
    events = map(Event, posts)
    return events

def print_events(events):
    "Prints the XML that corresponds to the events"
    print Event.xml_head
    print Event.start_tag
    events = filter(lambda e: e.allow, events)
    for event in events:
        print event.xml()
    print Event.end_tag


def type_by_date(events, start="Jan 1, 2000", end="Jan 1, 2020", type="something", days=10):
    "Activates events by date"
    start = datetime.strptime(start, "%b %d, %Y")
    end = datetime.strptime(end, "%b %d, %Y",)
    assert end > start, "%s, %s" % (start, end)
    for event in events:
        if start < event.date < end and event.type != SKIP:
            event.type = type

def keep(events):
    for event in events:
        event.allow = event.type

def filter_events(events):
    "This filters events by various parameters"


    # default view is that of a new user
    type_by_date(events, start="Aug 1, 2000", end="Dec 1, 2020", type="newuser")

    # when to show activity by all users
    type_by_date(events, start="Jan 1, 2000", end="May 1, 2010", type="users")
    type_by_date(events, start="Jan 1, 2011", end="Feb 1, 2011", type="users")
    type_by_date(events, start="Jan 1, 2012", end="Feb 1, 2012", type="users")
    type_by_date(events, start="Jan 1, 2013", end="Feb 1, 2013", type="users")


    # when to show activity by moderators
    type_by_date(events, start="May 1, 2011", end="Jun 1, 2011", type="moderator")
    type_by_date(events, start="May 1, 2012", end="Jun 1, 2012", type="moderator")
    type_by_date(events, start="May 1, 2013", end="Jun 1, 2013", type="moderator")


    # skip to action
    #type_by_date(events, end="Aug 1, 2010", type=SKIP)

    # moderator transitions
    #type_by_date(events, start="Apr 1, 2010", end="May 1, 2010", type=SKIP)

    # user transitions
    #type_by_date(events, start="Jun 1, 2010", end="Jul 1, 2010", type=SKIP)

    # remove the rest
    #type_by_date(events, start="Aug 1, 2010", type=SKIP)

    #type_by_date(events, start="Apr 10, 2010", end="Apr 20, 2010", type="skip")


    for event in events:

        if event.type == "moderator":
            event.allow = event.author.profile.can_moderate

        if event.type == "newuser":
            elapsed = event.date - event.author.date_joined
            event.allow = elapsed < timedelta(days=30)

        if event.type == "skip":
            event.allow = False

    #type_by_date(events, start="Jan 1, 2010", end="Feb 10, 2010", type="newusers")

    #keep(events)

    return events

def generate():
    LIMIT=None
    events = fetch(limit=LIMIT)
    events = filter_events(events)
    print_events(events)

if __name__ == "__main__":
    generate()