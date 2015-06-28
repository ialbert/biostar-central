import logging

from django.test import TestCase

from faker import Factory

from biostar3.forum import auth
from biostar3.forum.models import User, Post

logging.disable(logging.INFO)

f = Factory.create()

def create_post():
    author = User.objects.create(name=f.name(), email=f.email())
    auth.add_user_attributes(author)

    post_data = dict(
        title = f.sentence(),
        content = f.text(),
        tags = ""
    )
    post = auth.create_toplevel_post(user=author, data=post_data)

    return author, post

class ModerationTests(TestCase):
    def test_user_cant_restore(self):
        "A normal user can't restore a closed post"

        normal_user = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(normal_user)

        author, post = create_post()
        post.status = Post.CLOSED
        post.save()

        self.assertFalse("restore" in post.moderation_actions(normal_user))

    def test_author_cant_restore(self):
        "The author can't restore a closed post"

        author, post = create_post()
        post.status = Post.CLOSED
        post.save()

        self.assertFalse("restore" in post.moderation_actions(author))

    def test_moderator_can_restore(self):
        "A moderator can restore a closed post"

        moderator_user = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator_user)

        author, post = create_post()
        post.status = Post.CLOSED
        post.save()

        self.assertTrue("restore" in post.moderation_actions(moderator_user))

    def test_admin_can_restore(self):
        "An admin can restore a closed post"

        admin_user = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin_user)

        author, post = create_post()
        post.status = Post.CLOSED
        post.save()

        self.assertTrue("restore" in post.moderation_actions(admin_user))

    def test_can_close_top_level_post(self):
        "Can close a top-level post"

        moderator_user = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator_user)

        author, post = create_post()

        self.assertTrue("close" in post.moderation_actions(moderator_user))

    def test_cant_close_child_post(self):
        "Can't close a child post"

        moderator_user = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator_user)

        author, post = create_post()
        reply = auth.create_content_post(user=author, parent=post, content=f.sentence())

        self.assertFalse("close" in reply.moderation_actions(moderator_user))

    def test_can_mark_duplicate_top_level_post(self):
        "Can mark a top-level post as a duplicate"

        moderator_user = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator_user)

        author, post = create_post()

        self.assertTrue("duplicate" in post.moderation_actions(moderator_user))

    def test_cant_mark_duplicate_child_post(self):
        "Can't mark child posts as a duplicate"

        moderator_user = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator_user)

        author, post = create_post()
        reply = auth.create_content_post(user=author, parent=post, content=f.sentence())

        self.assertFalse("duplicate" in reply.moderation_actions(moderator_user))

    def test_cant_move_top_level_post(self):
        "Can't move top-level posts to other post types"

        moderator_user = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator_user)

        author, post = create_post()

        self.assertFalse("move_to_answer" in post.moderation_actions(moderator_user))
        self.assertFalse("move_to_comment" in post.moderation_actions(moderator_user))
        self.assertFalse("move_to_question" in post.moderation_actions(moderator_user))

    def test_can_move_child_post(self):
        "Can move child posts to other post types"

        moderator_user = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator_user)

        author, post = create_post()
        reply = auth.create_content_post(user=author, parent=post, content=f.sentence())

        self.assertTrue("move_to_answer" in reply.moderation_actions(moderator_user))
        self.assertTrue("move_to_comment" in reply.moderation_actions(moderator_user))
        self.assertTrue("move_to_question" in reply.moderation_actions(moderator_user))
