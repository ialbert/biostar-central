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
from biostar.apps.planet.models import BlogPost

from collections import defaultdict
from biostar.awards import create_user_award, check_user_profile

logger = logging.getLogger(__name__)

from allauth.socialaccount.adapter import DefaultSocialAccountAdapter


def get_ip(request):
    ip1 = request.META.get('REMOTE_ADDR', '')
    ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    ip = ip1 or ip2 or '0.0.0.0'
    return ip


class AutoSignupAdapter(DefaultSocialAccountAdapter):
    def pre_social_login(self, request, sociallogin):

        # This social login already exists.
        if sociallogin.is_existing:
            return

        # Check if we could/should connect it.
        try:
            email = sociallogin.account.extra_data.get('email')
            #verified = sociallogin.account.extra_data.get('verified_email')
            if email:
                user = User.objects.get(email=email)
                sociallogin.connect(request, user)
        except User.DoesNotExist:
            pass


class ExternalAuth(object):
    '''
    This is an "autentication" that relies on the user being valid.
    We're just following the Django interfaces here.
    '''

    def authenticate(self, email, valid=False):
        # Check the username/password and return a User.
        if valid:
            user = User.objects.get(email=email)
            user.backend = "%s.%s" % (__name__, self.__class__.__name__)
            print(user.backend)
            return user
        else:
            return None

    def get_user(self, user_id):
        try:
            return User.objects.get(pk=user_id)
        except User.DoesNotExist:
            return None

def valid_external_login(request):
    "Attempts to perform an external login"

    for name, key in settings.EXTERNAL_AUTH:
        value = request.COOKIES.get(name)
        if value:
            try:
                email, digest1 = value.split(":")
                digest2 = hmac.new(key, email).hexdigest()
                valid = (digest1 == digest2)
                if not valid:
                    raise Exception("digests do not match")
            except Exception as exc:
                logger.error(exc)
                return False

            # If we made it this far the data is valid.
            user, flag = User.objects.get_or_create(email=email)
            if flag:
                logger.info("created user %s" % user.email)

            # Authenticate with local info.
            user = ExternalAuth().authenticate(email=user.email, valid=valid)
            login(request=request, user=user)
            return True

    return False


SESSION_KEY, ANON_USER = settings.SESSION_KEY, "anon-user"


def get_counts(request, weeks=settings.COUNT_INTERVAL_WEEKS):
    "Returns the number of counts for each post type in the interval that has passed"
    user = request.user
    now = const.now()

    # Authenticated users get counts since their last login.
    if user.is_authenticated():
        since = user.profile.last_login
    else:
        since = now - timedelta(weeks=weeks)

    # This fetches the posts since last login.
    posts = Post.objects.filter(type__in=Post.TOP_LEVEL, status=Post.OPEN, creation_date__gt=since).order_by(
        '-id').only("id").prefetch_related("tag_set")
    posts = posts[:200]
    counts = defaultdict(int)

    # How many news posts.
    counts['latest'] = len(posts)

    # Produce counts per tag.
    for post in posts:
        for tag in post.tag_set.all():
            counts[tag.name] += 1

    # Fill in the unanswered counts.
    counts['open'] = Post.objects.filter(type=Post.QUESTION, reply_count=0, status=Post.OPEN,
                                         creation_date__gt=since).count()

    # How many new planet posts
    counts['planet'] = BlogPost.objects.filter(insert_date__gt=since).count()

    # Compute a few more counts for the user.
    if user.is_authenticated():
        # These are the new messages since the last login.
        counts['messages'] = Message.objects.filter(user=user, unread=True, sent_at__gt=since).count()

        # These are the new votes since the last login.
        counts['votes'] = Vote.objects.filter(post__author=user, date__gt=since).count()

    return counts


class Visit(object):
    """
    Sets visit specific parameters on objects.
    """

    def process_request(self, request, weeks=settings.COUNT_INTERVAL_WEEKS):
        global SESSION_KEY, ANON_USER

        user, session = request.user, request.session

        # Suspended users are logged out immediately.
        if user.is_authenticated() and user.is_suspended:
            logout(request)
            messages.error(request, 'Sorry, this account has been suspended. Please contact the administrators.')

        # Add attributes to anonymous users.
        if not user.is_authenticated():

            # This attribute is required inside templates.
            user.is_moderator = user.is_admin = False

            # Check external logins.
            if settings.EXTERNAL_AUTH and valid_external_login(request):
                messages.success(request, "Login completed")

            # We do this to detect when an anonymous session turns into a logged in one.
            if ANON_USER not in session:
                session[ANON_USER] = True

        # User attributes that refresh at given intervals.
        if user.is_authenticated():

            # The time between two count refreshes.
            elapsed = (const.now() - user.profile.last_login).seconds

            # The user has an anonymous session already.
            # Update the user login data now.
            if ANON_USER in session:
                del session[ANON_USER]
                elapsed = settings.SESSION_UPDATE_SECONDS + 1

            # The user session will be updated.
            if elapsed > settings.SESSION_UPDATE_SECONDS:
                # Set the last login time.
                Profile.objects.filter(user_id=user.id).update(last_login=const.now())

                # Compute the counts.
                counts = get_counts(request)

                # Store the counts in the session for later use.
                session[SESSION_KEY] = counts

                # Create user awards if possible.
                create_user_award.delay(user=user)

                # check user and fill in details
                check_user_profile.delay(ip=get_ip(request), user=user)


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
