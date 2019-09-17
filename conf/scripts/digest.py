import argparse
from datetime import timedelta
from django.template import loader
from biostar.forum.models import Post
from biostar.accounts.models import Profile
from biostar.emailer.tasks import send_email
from biostar.accounts import util


def send_digests(days=1, subject=""):

    prefs_map = {1: Profile.DAILY_DIGEST, 7: Profile.WEEKLY_DIGEST, 30: Profile.MONTHLY_DIGEST}
    # Send digests for posts within the last x days.
    delta = util.now() - timedelta(days=days)

    # Get users with the appropriate digest prefs.
    digest_prefs = prefs_map.get(days, Profile.DAILY_DIGEST)
    users = Profile.objects.filter(digest_prefs=digest_prefs)

    # Filter for users whose last digest sent is greater than or equal to the amount of days given.
    users = users.filter(last_digest__gte=delta)

    # Fetch posts within the last x amount of days
    posts = Post.objects.filter(lastedit_date__gt=delta)

    email_template = loader.get_template("messages/digest.html")
    context = dict(subject=subject, posts=posts)

    # Queue and send digest emails.
    emails = users.values_list('email', flat=True)
    send_email(template_name=email_template, extra_context=context, recipient_list=emails)

    # Update the last edit date to today
    users.update(last_digest=util.now())

    return


def test():
    # Testing the cron tasks
    print("in main")

    print("TESTING")

    import os
    test = os.path.abspath(os.path.join(os.path.dirname(__file__), 'tmp'))
    print(test)

    open(test, 'w').write('TESTING this cron')

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--daily', dest='daily', action='store_true', help='Send daily digests.')
    parser.add_argument('--weekly', dest='weekly', action='store_true', help='Send weekly digests.')
    parser.add_argument('--monthly', dest='monthly', action='store_true', help='Send monthly digests.')

    test()
    pass

    args = parser.parse_args()

    if args.daily:
        send_digests(days=1, subject="Daily digest")
    elif args.weekly:
        send_digests(days=7, subject="Weekly digest")
    elif args.monthly:
        send_digests(days=30, subject="Monthly digest")
