"""
Notification related tests.

These will execute when you run "manage.py test".
"""

import logging
from django.conf import settings
from biostar.apps.users.models import User, Profile
from biostar.apps.posts.models import Post, Subscription
from biostar.apps.messages.models import Message
from django.core import mail

from django.test import TestCase

logging.disable(logging.CRITICAL)


class NoteTest(TestCase):

    def test_send_email(self):
        "Testing email sending"
        mail.send_mail('Subject here', 'Here is the message.',
            'from@example.com', ['to@example.com'],
            fail_silently=False)
        self.assertEquals(len(mail.outbox), 1)
        self.assertEquals(mail.outbox[0].subject, 'Subject here')

    def test_note_creation(self):
        "Testing notifications"

        eq = self.assertEqual

        # Create some users
        title = "Test"
        emails = ["john@this.edu", "jane@this.edu", "bob@this.edu", "alice@this.edu",
                  "bill@this.edu", "jeff@this.edu" ]

        EMAIL_COUNT = len(emails)

        users, posts, user = [], [], None
        parent = None
        for email in emails:

            # Create users.
            user = User.objects.create(email=email)
            users.append(user)

        # Create a question.
        first = users[0]
        post = Post(title=title, author=first, type=Post.QUESTION)
        post.save()

        answers = []
        for user in users:
            # Every user adds an answer.
            answ = Post(author=user, type=Post.ANSWER, parent=post)
            answ.save()
            answers.append(answ)

        # A default admin user is added.
        eq(EMAIL_COUNT + 1, User.objects.all().count())

        # Total number of posts
        eq(EMAIL_COUNT + 1, Post.objects.all().count())

        # Every user has one subscription to the main post
        eq(EMAIL_COUNT, Subscription.objects.all().count())

        # Each user has a messages for content posted after
        # they started following the thread.
        for index, user in enumerate(users):
            mesg_c = Message.objects.filter(user=user).count()
            eq (mesg_c, EMAIL_COUNT - index - 1)

        '''
        COMMENT_COUNT = 3
        comments = []
        for i in range(COMMENT_COUNT):
            com = Post(author=last, type=Post.COMMENT, parent=answers[0])
            com.save()
            comments.append(com)

        # This is the total number of posts created.
        eq(EMAIL_COUNT + COMMENT_COUNT + ANSWER_COUNT, Post.objects.all().count())

        ansc = Post.objects.filter(type=Post.ANSWER).count()
        eq(ansc, 1)


        comc = Post.objects.filter(type=Post.COMMENT).count()
        eq(comc, 3)

        # Does not matter which comment we pick, each person should
        # have a subscription to the root.
        for sub in Subscription.objects.all():
            print sub.user, sub.post.root

        subsc = Subscription.objects.get_subs(post=comments[-1]).count()
        eq(subsc, EMAIL_COUNT)

        # Check how many messages do users have.
        # A user gets messages for posts they did not create.


        # now test the number of messages that they have
        for user in users:
            mesg_c = Message.objects.filter(user=user).count()
            post_c = Post.objects.exclude(author=user).count()

            print (user.email, post_c, mesg_c)


            #eq(num,
        #
        #    num = Message.objects.filter(user__email=email).count()
        #    print (email, index, num)
        #    eq(num, email_count + COMMENTS - index)

        '''

