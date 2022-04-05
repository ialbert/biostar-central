import logging
import inspect
import functools
from urllib.request import urlopen, Request
from functools import partial
import json
import mistune
from django.conf import settings
from django.template import loader
from biostar.utils.decorators import task


#
# Do not use logging in tasks! Deadlocking may occur!
#
# https://github.com/unbit/uwsgi/issues/1369
#

def message(msg, level=0):
    print(f"{msg}")


@task
def detect_location(ip, user_id):
    """
    Fills the user location based on url.
    """
    from biostar.accounts.models import Profile

    msg = f"location check for \tid={user_id}\tip={ip}"
    # The lookup needs to be turned on.
    if not settings.LOCATION_LOOKUP:
        message(f"skip {msg}")
        return

    message(f"{msg}")

    # Get the profile for the user
    profile = Profile.objects.filter(user__id=user_id).first()
    template = "messages/location-set.md"

    # Skip value if it has the word unknown in it
    def get(data, attr):
        value = data.get(attr, '')
        return "" if "unknown" in value.lower() else value.title()

    # Check and log location.
    if not profile.location:
        try:
            url = f"http://ip-api.com/json/{ip}"
            message(url)
            message(f"{ip}, {profile.user}, {url}")
            req = Request(url=url, headers={'User-Agent': 'Mozilla/5.0'})
            resp = urlopen(req, timeout=3).read()
            data = json.loads(resp)

            # Add the userid to the json
            data['uid'] = user_id

            # Log the return data.
            message(data)

            # city = get(data, "city")
            # region = get(data, "regionName")
            # country = get(data, "country")
            location = get(data, "country")

            msg = f"location result for \tid={user_id}\tip={ip}\tloc={location}"
            if location:
                Profile.objects.filter(user__id=user_id).update(location=location)
                context = dict(profile=profile, location=location)
                create_messages(template=template, user_ids=[user_id], extra_context=context)
                message(f"updated profile: {msg}")
            else:
                message(f"empty location: {msg}")

        except Exception as exc:

            message(exc)


@task
def verification_email(user_id):
    from biostar.accounts import auth, models

    user = models.User.objects.filter(id=user_id).first()

    auth.send_verification_email(user=user)
    return


@task
def create_messages(template, user_ids, sender=None, extra_context={}):
    """
    Create batch message from sender to a given recipient_list
    """
    from biostar.accounts.models import User, Message, MessageBody

    rec_list = User.objects.filter(id__in=user_ids)
    # Get the sender
    name, email = settings.ADMINS[0]
    sender = sender or User.objects.filter(email=email).first() or User.objects.filter(is_superuser=True).first()
    # Load the template and context
    tmpl = loader.get_template(template_name=template)

    # Default context added to each template.
    context = dict(sender=sender, domain=settings.SITE_DOMAIN, protocol=settings.PROTOCOL)

    context.update(extra_context)
    body = tmpl.render(context)
    html = mistune.markdown(body, escape=False)
    body = MessageBody.objects.create(body=body, html=html)

    for rec in rec_list:
        Message.objects.create(sender=sender, recipient=rec, body=body)
