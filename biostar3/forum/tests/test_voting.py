import logging

from django.test import TestCase

from faker import Factory

from biostar3.forum import auth, ajax
from biostar3.forum.models import User, Post, Vote

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

class VotingTests(TestCase):
    def test_author_does_not_gain_score(self):
        "The post author does not gain score by upvoting their own post"

        author, post = create_post()
        original_score = author.score

        ajax.perform_vote(post, author, Vote.UP)

        # Reload author from database
        author = User.objects.get(id=author.id)

        self.assertEqual(author.score, original_score)

    def test_vote_increases_score(self):
        "Upvoting causes the author's score to increase"

        author, post = create_post()
        original_score = author.score
        other_user = User.objects.create(name=f.name(), email=f.email())

        ajax.perform_vote(post, other_user, Vote.UP)

        # Reload author from database
        author = User.objects.get(id=author.id)

        self.assertGreater(author.score, original_score)

    def test_unvote_decreases_score(self):
        "Un-upvoting causes the author's score to decrease"

        author, post = create_post()
        original_score = author.score
        other_user = User.objects.create(name=f.name(), email=f.email())

        ajax.perform_vote(post, other_user, Vote.UP)

        # Reload author from database
        author = User.objects.get(id=author.id)
        upvoted_score = author.score

        ajax.perform_vote(post, other_user, Vote.UP)

        # Reload author from database
        author = User.objects.get(id=author.id)

        self.assertLess(author.score, upvoted_score)
        self.assertEqual(author.score, original_score)

    def test_upvoting_child_increases_thread_score(self):
        "Upvoting a child post causes the root's thread score to increase"

        author, root_post = create_post()
        reply = auth.create_content_post(user=author, parent=root_post, content=f.sentence())
        original_thread_score = root_post.thread_score
        other_user = User.objects.create(name=f.name(), email=f.email())

        ajax.perform_vote(reply, other_user, Vote.UP)

        # Reload root_post from database
        root_post = Post.objects.get(id=root_post.id)

        self.assertEqual(root_post.thread_score, original_thread_score + 1)

    def test_accepting_sets_has_accepted(self):
        "Accepting an answer sets has_accepted on the answer and question"

        asker, question = create_post()
        answerer = User.objects.create(name=f.name(), email=f.email())
        answer = auth.create_content_post(user=answerer, parent=question, content=f.sentence())

        self.assertFalse(question.has_accepted)
        self.assertFalse(answer.has_accepted)

        ajax.perform_vote(answer, asker, Vote.ACCEPT)

        # Reload posts from database
        question = Post.objects.get(id=question.id)
        answer = Post.objects.get(id=answer.id)

        self.assertTrue(question.has_accepted)
        self.assertTrue(answer.has_accepted)

    def test_multiple_accept(self):
        "Accepting multiple answers is possible"

        asker, question = create_post()
        answerer = User.objects.create(name=f.name(), email=f.email())
        answer1 = auth.create_content_post(user=answerer, parent=question, content=f.sentence())
        answer2 = auth.create_content_post(user=answerer, parent=question, content=f.sentence())

        ajax.perform_vote(answer1, asker, Vote.ACCEPT)

        # Reload posts from database
        question = Post.objects.get(id=question.id)
        answer1 = Post.objects.get(id=answer1.id)
        answer2 = Post.objects.get(id=answer2.id)

        self.assertTrue(question.has_accepted)
        self.assertTrue(answer1.has_accepted)
        self.assertFalse(answer2.has_accepted)

        ajax.perform_vote(answer2, asker, Vote.ACCEPT)

        # Reload posts from database
        question = Post.objects.get(id=question.id)
        answer1 = Post.objects.get(id=answer1.id)
        answer2 = Post.objects.get(id=answer2.id)

        self.assertTrue(question.has_accepted)
        self.assertTrue(answer1.has_accepted)
        self.assertTrue(answer2.has_accepted)
