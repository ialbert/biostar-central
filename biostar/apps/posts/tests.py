"""
Post related tests.

These will execute when you run "manage.py test".
"""
from __future__ import print_function, unicode_literals, absolute_import, division

import logging
from django.conf import settings
from biostar.apps.users.models import User, Profile
from biostar.apps.posts.models import Post, Subscription, Tag
from biostar.apps.messages.models import Message

from django.test import TestCase

logging.disable(logging.INFO)

TEST_POST = '''
# Hello world


Linkify urls in the test post http://www.biostars.org

1. One
2. Two
3. Three

Do not linkify posts in code:

    ls *
    cat data.fq | grep "HWUS" | wc -l
    curl http://www.biostars.org

Reformat links to internal content:

* Post: http://www.lvh.me:8080/p/22/
* User: http://www.lvh.me:8080/u/56/

Embed tweets:

https://twitter.com/Linux/status/2311234267

Embed youtube:

https://www.youtube.com/watch?v=kfvxmEuC7bU

Embed gist:

https://gist.github.com/ialbert/ae46c5f51d63cdf2d0d2

'''

class PostTest(TestCase):

    def test_tagging(self):
        "Testing tagging."
        eq = self.assertEqual

        eq(0, Tag.objects.all().count() )

        # Create an admin user and a post.
        title = "Hello Posts!"
        email = "john@this.edu"
        jane = User.objects.create(email=email)
        html = "<b>Hello World!</b>"
        post = Post(title=title, author=jane, type=Post.FORUM, content=html)
        post.save()
        post.add_tags("t1,t2, t3")

        eq(3, Tag.objects.all().count())

        post = Post(title=title, author=jane, type=Post.FORUM, content=html)
        post.save()
        post.add_tags("t1, t2, t3, t2, t1, t1")

        t1 = Tag.objects.get(name="t1")
        t3 = Tag.objects.get(name="t3")

        eq(2, t1.count)
        eq(2, t3.count)

        post.add_tags("t2 t4")

        t1 = Tag.objects.get(name="t1")
        t3 = Tag.objects.get(name="t3")

        eq(1, t1.count)
        eq(1, t3.count)

    def test_post_creation(self):
        "Testing post creation."
        eq = self.assertEqual

        # Create an admin user and a post.
        title = "Hello Posts!"
        email = "john@this.edu"
        jane = User.objects.create(email=email)
        html = "<b>Hello World!</b>"
        post = Post(title=title, author=jane, type=Post.FORUM, content=html)
        post.save()

        # Get the object fresh.
        post = Post.objects.get(pk=post.id)

        eq(post.type, Post.FORUM)
        eq(post.root, post)
        eq(post.parent, post)

        # Subscriptions are automatically created
        sub = Subscription.objects.get(user=jane)
        eq(sub.user, jane)
        eq(sub.post, post)

        # A new post triggers a message to the author.
        email = "jane@this.edu"
        john = User.objects.create(email=email)
        answer = Post(author=john, parent=post, type=Post.ANSWER)
        answer.save()

        eq(answer.root, post)
        eq(answer.parent, post)
        eq(answer.type, Post.ANSWER)

        # Add comment. The parent will override the post type.
        email = "bob@this.edu"
        bob = User.objects.create(email=email)
        comment = Post(author=bob, type=Post.FORUM, parent=answer)
        comment.save()

        eq(comment.root, post)
        eq(comment.parent, answer)
        eq(comment.type, Post.COMMENT)

        # Everyone posting in a thread gets a subscription to the root post of the
        subs = Subscription.objects.filter(post=post)
        eq(len(subs), 3)

TEST_CONTENT_EMBEDDING ="""
<p>Gist links may be formatted</p>

<pre>
https://gist.github.com/ialbert/ae46c5f51d63cdf2d0d2</pre>

<p>or embedded:</p>

<p>https://gist.github.com/ialbert/ae46c5f51d63cdf2d0d2</p>

<p>Video links may be formatted</p>

<pre>
http://www.youtube.com/watch?v=_cDaX0xJPvI</pre>

<p>or embedded:</p>

<p>http://www.youtube.com/watch?v=_cDaX0xJPvI</p>

<p>Internal links are recognized:</p>

<pre>
http://test.biostars.org/u/2/</pre>


<p>vs&nbsp;http://test.biostars.org/u/2/</p>
<p>Similarly&nbsp;</p>

<pre>
http://test.biostars.org/p/2/</pre>

<p>versus&nbsp;http://test.biostars.org/p/2/</p>

<p>&nbsp;</p>
"""