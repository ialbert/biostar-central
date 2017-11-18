from django.contrib import admin
from .models import Job, Analysis




class JobAdmin(admin.ModelAdmin):

    search_fields = ('name', 'owner__first_name', 'owner__email', 'state')
    list_display = ("name", "state","start_date", "security","date")
    list_filter = ("state", "security")


class AnalysisAdmin(admin.ModelAdmin):

    search_fields = ('name', 'text', 'owner__first_name', 'owner__email' )
    list_display = ("name", "date", "auth")
    list_filter = ( "state", "auth")


admin.site.register(Job, JobAdmin)
admin.site.register(Analysis, AnalysisAdmin)