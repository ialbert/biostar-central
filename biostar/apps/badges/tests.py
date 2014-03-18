from django.test import TestCase
from .models import Award, Badge
from biostar.apps.users.models import User

# Create your tests here.

class AwardTest(TestCase):
    email = "janedoe@site.com"

    def setUp(self):
        from biostar import awards
        awards.init_awards()
        User.objects.create(email=self.email)

    def test_user_badge(self):
        from biostar import awards
        eq = self.assertEqual

        award_count = lambda: Award.objects.all().count()

        for a in Award.objects.all():
            print a

        eq(0, award_count())

        jane = User.objects.get(email=self.email)
        awards.create_user_award(jane)

        # No award for new user.
        eq(0, award_count())

        jane.profile.info = "A" * 1000
        jane.save()

        # Check for the autobiographer award.
        awards.create_user_award(jane)

        eq(1, award_count())
