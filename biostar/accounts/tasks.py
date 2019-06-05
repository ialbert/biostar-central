import logging
import hjson

from urllib.request import urlopen, Request
from django.conf import settings
from django.template import loader

from biostar.accounts import models

logger = logging.getLogger('biostar')


def detect_location(request, user):
    """
    Fills the user location based on url.
    """

    # Get the ip information
    ip1 = request.META.get('REMOTE_ADDR', '')
    ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    ip = ip1 or ip2 or '0.0.0.0'

    logger.info(f"location-check\tid={user.id}\tip={ip}")

    # Don't hammer the servers when testing
    if settings.DEBUG:
        return

    # Check and log location.
    if user.is_authenticated and not user.profile.location:

        try:
            url = f"http://api.hostip.info/get_json.php?ip={ip}"
            logger.debug(f"{ip}, {user}, {url}")
            req = Request(url=url, headers={'User-Agent': 'Mozilla/5.0'})
            resp = urlopen(req, timeout=3).read()
            data = hjson.loads(resp)
            city = data.get("city", '')
            country = data.get('country_name', '').title()
            location = city or country
            location = "localhost" if ip in ('127.0.0.1') else location
            if "unknown" not in location.lower():
                models.Profile.objects.filter(user=user).update(location=location)
                logger.info(f"location-set\tid={user.id}\tip={ip}\tloc={location}")

        except Exception as exc:
            logger.error(exc)


def create_messages(template, rec_list, sender=None, extra_context={}, subject=""):

    """
    Create batch message from sender to a given recipient_list
    """
    # Get the sender
    name, email = settings.ADMINS[0]
    sender = sender or models.User.objects.filter(email=email).first()

    # Load the template and context
    tmpl = loader.get_template(template_name=template)
    context = dict(sender=sender, subject=subject)
    context.update(extra_context)

    body = tmpl.render(context)

    msgs = []
    for rec in rec_list:
        msg = models.Message.objects.create(sender=sender, recipient=rec, subject=subject, body=body)
        msgs.append(msg)

    return msgs
