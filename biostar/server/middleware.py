__author__ = 'ialbert'
from django.contrib import messages
from django.conf import settings
import hmac, logging, re
from datetime import timedelta
from django.contrib.auth import authenticate, login, logout
from biostar.apps.users.models import User, Profile
from biostar import const
from django.core.cache import cache
from biostar.apps.posts.models import Post, Vote
from biostar.apps.messages.models import Message

from collections import defaultdict

logger = logging.getLogger(__name__)

def valid_external_login(request):
    "Attempts to perform an external login"

    for name, key in settings.EXTERNAL_AUTH:
        value = request.COOKIES.get(name)
        if value:
            try:
                email, digest1 = value.split(":")
                digest2 = hmac.new(key, email).hexdigest()
                if digest1 != digest2:
                    raise Exception("digests do not match: %s vs %s" % (digest1, digest2))
            except Exception, exc:
                logger.error(exc)
                return None

            # If we made it this far the data is valid.
            password = settings.SECRET_KEY + email

            user, flag = User.objects.get_or_create(email=email)
            if flag:
                logger.info("created user %s" % user.email)

            # Regular login is not allowed.
            user.set_password(password)
            user.save()

            # Authenticate.
            user = authenticate(username=user.email, password=password)
            login(request=request, user=user)
            return True

    return False

SESSION_KEY = settings.SESSION_KEY

def get_counts(request, weeks=settings.COUNT_INTERVAL_WEEKS):
    "Returns the number of counts for each post type in the interval that has passed"
    user = request.user
    now = const.now()

    if user.is_authenticated():
        since = user.profile.last_login
    else:
        since = now - timedelta(weeks=weeks)

    posts = Post.objects.filter(type__in=Post.TOP_LEVEL, status=Post.OPEN, creation_date__gt=since).order_by('-id').only("id").prefetch_related("tag_set")
    counts = defaultdict(int)

    # How many news posts.
    counts['latest'] = len(posts)

    # Produce counts per tag.
    for post in posts:
        for tag in post.tag_set.all():
            counts[tag.name] += 1

    # Fill in the unanswered counts.
    counts['unanswered'] = Post.objects.filter(type=Post.QUESTION, reply_count=0, status=Post.OPEN, creation_date__gt=since).count()

    # Compute a few more counts for the user.
    if user.is_authenticated():
        # These are the new messages since the last login.
        counts['messages'] = Message.objects.filter(user=user, unread=None, sent_at__gt=since).count()

        # These are the new votes since the last login.
        counts['votes'] = Vote.objects.filter(post__author=user, date__gt=since).count()

    return counts


class Visit(object):
    """
    Sets visit specific parameters on objects.
    """

    def process_request(self, request, weeks=settings.COUNT_INTERVAL_WEEKS):
        global SESSION_KEY

        user, session = request.user, request.session

        # Suspended users are logged out immediately.
        if user.is_authenticated() and user.is_suspended:
            logout(request)
            messages.error(request, 'Sorry, this account has been suspended. Please contact the administrators.')
            return

        # Add attributes to anonymous users.
        if not user.is_authenticated():

            # This attribute is required inside templates.
            user.is_moderator = user.is_admin = False

            # Check external logins.
            if settings.EXTERNAL_AUTH and valid_external_login(request):
                messages.success(request, "Login completed")

        # User attributes that refresh at given intervals.
        if user.is_authenticated():

            # The time between two count refreshes.
            elapsed = (const.now() - user.profile.last_login).seconds

            if elapsed > settings.SESSION_UPDATE_SECONDS:
                # Set the last login time.
                Profile.objects.filter(user_id=user.id).update(last_login=const.now())

                # Compute the counts.
                counts = get_counts(request)

                # Store the counts in the session for later use.
                session[SESSION_KEY] = counts

        # Get the counts from the session or the cache.
        counts = session.get(SESSION_KEY) or cache.get(SESSION_KEY)

        # No sessions found, set the them into the session.
        if not counts:
            # Compute the counts
            counts = get_counts(request)

            # Put them into the session.
            session[SESSION_KEY] = counts

            # Store them in the cache for the next anonymous user.
            cache.set(SESSION_KEY, counts, settings.SESSION_UPDATE_SECONDS)

        # The session key must be present in the session
        assert SESSION_KEY in session

