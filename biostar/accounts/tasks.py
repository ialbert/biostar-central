import logging
import hjson

from urllib.request import urlopen, Request
from django.conf import settings
from django.template import loader

from biostar.utils.decorators import spool

logger = logging.getLogger('biostar')

@spool(pass_arguments=True)
def detect_location(ip, user_id):
    """
    Fills the user location based on url.
    """
    from biostar.accounts.models import Profile

    logger.info(f"location-check\tid={user_id}\tip={ip}")

    # Don't hammer the servers when testing
    if settings.DEBUG:
        return

    # Get the profile for the user
    profile = Profile.objects.filter(user__id=user_id).first()

    # Check and log location.
    if not profile.location and settings.LOCATION_LOOKUP:
        try:
            url = f"http://api.hostip.info/get_json.php?ip={ip}"
            logger.debug(f"{ip}, {profile.user}, {url}")
            req = Request(url=url, headers={'User-Agent': 'Mozilla/5.0'})
            resp = urlopen(req, timeout=3).read()
            data = hjson.loads(resp)
            city = data.get("city", '')
            country = data.get('country_name', '').title()
            location = city or country
            location = "localhost" if ip in ('127.0.0.1') else location
            if "unknown" not in location.lower():
                Profile.objects.filter(user=profile.user).update(location=location)
                logger.info(f"location-set\tid={profile.user.id}\tip={ip}\tloc={location}")

        except Exception as exc:
            logger.error(exc)

@spool(pass_arguments=True)
def create_messages(template, rec_list, sender=None, extra_context={}, subject=""):
    """
    Create batch message from sender to a given recipient_list
    """
    from biostar.accounts.models import User, Message

    # Get the sender
    name, email = settings.ADMINS[0]
    sender = sender or User.objects.filter(email=email).first()

    # Load the template and context
    tmpl = loader.get_template(template_name=template)
    context = dict(sender=sender, subject=subject)
    context.update(extra_context)

    body = tmpl.render(context)

    msgs = []
    for rec in rec_list:

        msg = Message.objects.create(sender=sender, recipient=rec, subject=subject, body=body)
        msgs.append(msg)

    return msgs
