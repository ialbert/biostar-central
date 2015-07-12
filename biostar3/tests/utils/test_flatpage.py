from django.test import TestCase

from biostar3.utils.flatpage import parse_metadata, add_one
from biostar3.forum.models import User, Post, FlatPage
from biostar3.forum import auth

import logging
logging.disable(logging.INFO)

from faker import Factory
f = Factory.create()

class ParseMetadataTests(TestCase):
    def test_single(self):
        content = "{# asdf = fooasdfbar #}"
        meta = parse_metadata(content)

        self.assertEqual(meta, dict(asdf="fooasdfbar"))

    def test_value_has_spaces(self):
        content = "{# asdf = foo asdf bar #}"
        meta = parse_metadata(content)

        self.assertEqual(meta, dict(asdf="foo asdf bar"))

    def test_multi(self):
        content = ( "{# asdf = fooasdfbar #}\n"
                    "{# lkjh = barasdffoo #}" )
        meta = parse_metadata(content)

        self.assertEqual(meta, dict(asdf="fooasdfbar", lkjh="barasdffoo"))

    def test_with_other_content(self):
        content = ( "qwertyuiop\n"
                    "{# asdf = fooasdfbar #}\n"
                    "{# lkjh = barasdffoo #}\n"
                    "mnbvcxz" )
        meta = parse_metadata(content)

        self.assertEqual(meta, dict(asdf="fooasdfbar", lkjh="barasdffoo"))

    def test_ignore_empty_block(self):
        content = ( "{# #}" )
        meta = parse_metadata(content)

        self.assertEqual(meta, dict())

class AddOneTests(TestCase):
    def test_empty_should_error(self):
        content = ""
        user = User.objects.create(name=f.name(), email=f.email())
        path = "foo/bar.md"
        errors = add_one(user, content, path)

        self.assertIn("no content in {}".format(path), errors)

    def test_no_slug_should_error(self):
        content = "fooasdfbar"
        user = User.objects.create(name=f.name(), email=f.email())
        path = "foo/bar.md"
        errors = add_one(user, content, path)

        self.assertIn("slug field is missing from {}".format(path), errors)

    def test_no_title_should_error(self):
        content = "fooasdfbar"
        user = User.objects.create(name=f.name(), email=f.email())
        path = "foo/bar.md"
        errors = add_one(user, content, path)

        self.assertIn("title field is missing from {}".format(path), errors)

    def test_error_when_slug_exists_and_not_update(self):
        slug = "fooasdfbar"
        content = ( "{{# slug = {slug} #}}\n"
                    "{{# title = Fooasdfbar #}}\n"
                    "fooasdfbar" ).format(slug=slug)
        user = User.objects.create(name=f.name(), email=f.email())
        path = "foo/bar.md"

        old_post = auth.create_toplevel_post(data=dict(title=f.sentence(), type=Post.PAGE, content=f.sentence()), user=user)
        FlatPage.objects.create(post=old_post, slug=slug)

        errors = add_one(user, content, path)

        self.assertIn("slug {} already exists from {}".format(slug, path), errors)
        self.assertEqual(FlatPage.objects.filter(slug=slug).count(), 1)

    def test_does_create(self):
        slug = "fooasdfbar"
        content = ( "{{# slug = {slug} #}}\n"
                    "{{# title = Fooasdfbar #}}\n"
                    "fooasdfbar" ).format(slug=slug)
        user = User.objects.create(name=f.name(), email=f.email())
        path = "foo/bar.md"
        errors = add_one(user, content, path)

        self.assertFalse(errors)
        self.assertEqual(FlatPage.objects.filter(slug=slug).count(), 1)