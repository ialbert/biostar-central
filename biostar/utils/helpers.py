from django.contrib.messages.storage import fallback
from django.test import RequestFactory, client
from django.conf import settings
import logging
from datetime import datetime
from biostar import VERSION
import os
import uuid

logger = logging.getLogger('engine')


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def fake_request(url, data, user, method="POST", rmeta={}):
    "Make a fake request; defaults to POST."

    methods = {"POST": RequestFactory().post, "GET": RequestFactory().get,
               'PUT': RequestFactory().put}

    assert method in methods

    request = methods[method](url, data)

    # Mimic messaging system
    request.session = {}
    messages = fallback.FallbackStorage(request=request)
    request._messages = messages

    # Update the META info
    request.META.update(rmeta)

    request.user = user

    return request


def pg_dump(prog, pg_user, outdir, hourly=False,  **kwargs):
    # Get the full path to the directory.
    """
    """

    os.makedirs(outdir, exist_ok=True)

    pg_name = settings.DATABASE_NAME

    # Hourly database dumps have a simpler name so
    # that they overwrite each other.
    if hourly:
        # These names only include the hours.
        tstamp = datetime.now().strftime("hourly-%H")
    else:
        # These names include the date.
        tstamp = datetime.now().strftime("%Y-%m-%d")

    db_file = "%s-%s-%s.sql.gz" % (pg_name, VERSION, tstamp)
    db_file = os.path.abspath(os.path.join(outdir, db_file))

    params = dict(
        pg_user=pg_user,
        pg_name=pg_name,
        version=VERSION,
        db_file=db_file,
        prog=prog,
    )

    # logger.info("saving %(pg_name)s to %(db_file)s" % params)

    cmd = "%(prog)s -x -O -b -U %(pg_user)s %(pg_name)s | gzip > %(db_file)s" % params

    # Running the command
    logger.info("%s" % cmd)
    os.system(cmd)

    return


def get_ip(request):
    """
    Attempts to extract the IP number from the HTTP request headers.
    """
    # lower the
    key = settings.IP_HEADER_KEY
    meta = request.META

    # Lowercase keys
    lower = {k.lower(): v for k, v in request.META.items()}

    ip = meta.get(key, lower.get(key, '0.0.0.0'))

    return ip


def ip_triplet(request):
    """
    Attempt to extract first three number from ip adress.
    """
    oip = get_ip(request=request)
    ips = oip.split(".")[:-1]
    ip = ".".join(ips)
    return ip