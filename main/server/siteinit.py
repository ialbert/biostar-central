"""
Actions that need to take place when the site is initialized
"""
from main.server import models, const
from django.conf import settings

# upon startup create and set the password for admin if it does not exist
admin, created = models.User.objects.get_or_create(username='admin', email=settings.ADMINS[0][1], is_staff=True, is_superuser=True)
if created:
    admin.set_password(settings.SECRET_KEY)
    admin.save()
    print '*** created admin user'
    
admin = models.User.objects.get(username='admin')

# and a sanity check
assert admin.is_staff and admin.is_superuser

editors, created = models.Group.objects.get_or_create(name=const.EDITOR_GROUP)
if created:
    print '*** created group %s' % editors.name


