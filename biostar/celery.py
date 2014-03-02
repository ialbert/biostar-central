from __future__ import absolute_import
from django.conf import settings
from celery.utils.log import get_task_logger

from biostar import const
from datetime import timedelta

logger = get_task_logger(__name__)

from celery import Celery

app = Celery('biostar')

# Read the configuration from the config file.
app.config_from_object(settings.CELERY_CONFIG)

# Discover tasks in applications.
app.autodiscover_tasks(
    lambda: ["biostar.mailer"]
)


@app.task
def data_cleanup(days=1, weeks=20):
    from biostar.apps.posts.models import PostView
    from biostar.apps.messages.models import Message

    "Reduces messages and post views"
    past = const.now() - timedelta(days=days)
    query = PostView.objects.filter(date__lt=past)
    msg = "Deleting %s PostViews" % query.count()
    logger.info(msg)
    query.delete()

    # Reduce messages.
    since = const.now() - timedelta(weeks=weeks)
    query = Message.objects.filter(sent_at__lt=since)
    msg = "Deleting %s messages" % query.count()
    logger.info(msg)
    query.delete()


@app.task
def test(*args, **kwds):
    logger.info("*** executing task %s %s, %s" % (__name__, args, kwds))
    return 1000



