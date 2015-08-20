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

class PostModerationTests(TestCase):
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

class UserModerationTests(TestCase):
    def test_user_cant_moderate_users(self):
        "A normal user can't perform any user moderation actions"

        normal_user = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(normal_user)
        target_user = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(target_user)

        self.assertEqual(target_user.moderation_actions(normal_user), [])

    def test_moderator_cant_moderate_admin(self):
        "A moderator can't perform moderation actions on an admin"

        moderator = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator)
        admin = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin)

        self.assertEqual(admin.moderation_actions(moderator), [])

    def test_admin_cant_moderate_admin(self):
        "An admin can't perform moderation actions on an admin"

        admin1 = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin1)
        admin2 = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin2)

        self.assertEqual(admin2.moderation_actions(admin1), [])

    def test_staff_cant_moderate_admin(self):
        "Staff can't perform moderation actions on an admin"

        staff = User.objects.create(name=f.name(), email=f.email(), is_staff=True)
        auth.add_user_attributes(staff)
        admin = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin)

        self.assertEqual(admin.moderation_actions(staff), [])

    def test_moderator_cant_moderate_staff(self):
        "A moderator can't perform moderation actions on staff"

        moderator = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator)
        staff = User.objects.create(name=f.name(), email=f.email(), is_staff=True)
        auth.add_user_attributes(staff)

        self.assertEqual(staff.moderation_actions(moderator), [])

    def test_admin_cant_moderate_staff(self):
        "An admin can't perform moderation actions on staff"

        admin = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin)
        staff = User.objects.create(name=f.name(), email=f.email(), is_staff=True)
        auth.add_user_attributes(staff)

        self.assertEqual(staff.moderation_actions(admin), [])

    def test_staff_cant_moderate_staff(self):
        "Staff can't perform moderation actions on staff"

        staff1 = User.objects.create(name=f.name(), email=f.email(), is_staff=True)
        auth.add_user_attributes(staff1)
        staff2 = User.objects.create(name=f.name(), email=f.email(), is_staff=True)
        auth.add_user_attributes(staff2)

        self.assertEqual(staff2.moderation_actions(staff1), [])

    def test_moderator_can_suspend_user(self):
        "A moderator can suspend a user"

        moderator = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator)
        user = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user)

        self.assertFalse(user.is_suspended)

        self.assertTrue("suspend" in user.moderation_actions(moderator))

    def test_admin_can_suspend_user(self):
        "An admin can suspend a user"

        admin = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin)
        user = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user)

        self.assertFalse(user.is_suspended)

        self.assertTrue("suspend" in user.moderation_actions(admin))

    def test_staff_can_suspend_user(self):
        "Staff can suspend a user"

        staff = User.objects.create(name=f.name(), email=f.email(), is_staff=True)
        auth.add_user_attributes(staff)
        user = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user)

        self.assertFalse(user.is_suspended)

        self.assertTrue("suspend" in user.moderation_actions(staff))

    def test_moderator_cant_ban_user(self):
        "A moderator can't ban a user"

        moderator = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator)
        user = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user)

        self.assertFalse(user.is_suspended)

        self.assertFalse("ban" in user.moderation_actions(moderator))

    def test_admin_cant_ban_user(self):
        "An admin can't ban a user"

        admin = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin)
        user = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user)

        self.assertFalse(user.is_suspended)

        self.assertFalse("ban" in user.moderation_actions(admin))

    def test_staff_can_ban_user(self):
        "Staff can ban a user"

        staff = User.objects.create(name=f.name(), email=f.email(), is_staff=True)
        auth.add_user_attributes(staff)
        user = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user)

        self.assertFalse(user.is_suspended)

        self.assertTrue("ban" in user.moderation_actions(staff))

    def test_moderator_can_reinstate_suspended_user(self):
        "A moderator can reinstate a suspended user"

        moderator = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator)
        user = User.objects.create(name=f.name(), email=f.email(), status=User.SUSPENDED)
        auth.add_user_attributes(user)

        self.assertTrue(user.is_suspended)

        self.assertTrue("reinstate" in user.moderation_actions(moderator))

    def test_admin_can_reinstate_suspended_user(self):
        "An admin can reinstate a suspended user"

        admin = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin)
        user = User.objects.create(name=f.name(), email=f.email(), status=User.SUSPENDED)
        auth.add_user_attributes(user)

        self.assertTrue(user.is_suspended)

        self.assertTrue("reinstate" in user.moderation_actions(admin))

    def test_staff_can_reinstate_suspended_user(self):
        "Staff can reinstate a suspended user"

        staff = User.objects.create(name=f.name(), email=f.email(), is_staff=True)
        auth.add_user_attributes(staff)
        user = User.objects.create(name=f.name(), email=f.email(), status=User.SUSPENDED)
        auth.add_user_attributes(user)

        self.assertTrue(user.is_suspended)

        self.assertTrue("reinstate" in user.moderation_actions(staff))

    def test_moderator_cant_reinstate_banned_user(self):
        "A moderator can't reinstate a banned user"

        moderator = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator)
        user = User.objects.create(name=f.name(), email=f.email(), status=User.BANNED)
        auth.add_user_attributes(user)

        self.assertTrue(user.is_suspended)

        self.assertFalse("reinstate" in user.moderation_actions(moderator))

    def test_admin_cant_reinstate_banned_user(self):
        "An admin can't reinstate a banned user"

        admin = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin)
        user = User.objects.create(name=f.name(), email=f.email(), status=User.BANNED)
        auth.add_user_attributes(user)

        self.assertTrue(user.is_suspended)

        self.assertFalse("reinstate" in user.moderation_actions(admin))

    def test_staff_can_reinstate_banned_user(self):
        "Staff can reinstate a banned user"

        staff = User.objects.create(name=f.name(), email=f.email(), is_staff=True)
        auth.add_user_attributes(staff)
        user = User.objects.create(name=f.name(), email=f.email(), status=User.BANNED)
        auth.add_user_attributes(user)

        self.assertTrue(user.is_suspended)

        self.assertTrue("reinstate" in user.moderation_actions(staff))

    def test_cant_reinstate_normal_user(self):
        "Can't reinstate a user that isn't suspended or banned"

        moderator = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator)
        user = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user)

        self.assertFalse(user.is_suspended)

        self.assertFalse("reinstate" in user.moderation_actions(moderator))

    def test_moderator_cant_merge_users(self):
        "A moderator can't merge users"

        moderator = User.objects.create(name=f.name(), email=f.email(), type=User.MODERATOR)
        auth.add_user_attributes(moderator)
        user1 = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user1)
        user2 = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user2)

        self.assertFalse("merge" in user1.moderation_actions(moderator))
        self.assertFalse("merge" in user2.moderation_actions(moderator))

    def test_admin_cant_merge_users(self):
        "An admin can't merge users"

        admin = User.objects.create(name=f.name(), email=f.email(), type=User.ADMIN)
        auth.add_user_attributes(admin)
        user1 = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user1)
        user2 = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user2)

        self.assertFalse("merge" in user1.moderation_actions(admin))
        self.assertFalse("merge" in user2.moderation_actions(admin))

    def test_staff_can_merge_users(self):
        "Staff can merge users"

        staff = User.objects.create(name=f.name(), email=f.email(), is_staff=True)
        auth.add_user_attributes(staff)
        user1 = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user1)
        user2 = User.objects.create(name=f.name(), email=f.email())
        auth.add_user_attributes(user2)

        self.assertTrue("merge" in user1.moderation_actions(staff))
        self.assertTrue("merge" in user2.moderation_actions(staff))
