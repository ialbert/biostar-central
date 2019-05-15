from django.contrib import admin
from .models import EmailGroup,EmailAddress,Subscription

admin.site.register(EmailAddress)
admin.site.register(EmailGroup)
admin.site.register(Subscription)
