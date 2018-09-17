from django.contrib import admin
from .models import Job, Analysis, Project, Access
from django.forms import Textarea
from django.db.models import TextField


@admin.register(Access)
class AccessAdmin(admin.ModelAdmin):
    list_select_related = (
        'user',
        'project',
    )

    pass

@admin.register(Job)
class JobAdmin(admin.ModelAdmin):

    formfield_overrides = {
        TextField: {'widget': Textarea(attrs={'rows':20,'cols': 100})},
    }

    search_fields = ('name', 'owner__first_name', 'owner__email', 'state', "project__name",
                     "project__owner__first_name", "project__owner__email")
    list_display = ("name", "state","start_date", "security","date")
    list_filter = ("state", "security", "project__name", "deleted")


    fieldsets = (("Job Metadata",
                    {'fields': ("name","owner",'project',("uid","sticky"),
                                ("state", "security"), "image"),
                     "classes": ('extrapretty')}
                  ),

                 ("Optional Text Inputs",
                    {'fields':(("text", "html"),"summary"),
                     "classes": ("collapse",'extrapretty')}
                  ),

                 ("Run Time Settings",
                    {'fields':("json_text", "template"),
                     "classes": ("wide", 'extrapretty')},
                  ),
                 )


@admin.register(Analysis)
class AnalysisAdmin(admin.ModelAdmin):

    formfield_overrides = {
        TextField: {'widget': Textarea(attrs={'rows':20,'cols': 100})},
    }

    search_fields = ('name', 'text', 'owner__first_name', 'owner__email', "project__name",
                     "project__owner__first_name", "project__owner__email")
    list_display = ("name", "date", "security")
    list_filter = ("security", "project__name", "deleted")


    fieldsets = (("Analysis Metadata",
                    {'fields': ("name","owner",'project',("uid","sticky"),
                                ( "deleted", "security"), "image"),
                     "classes": ('extrapretty')}
                  ),

                 ("Optional Text Inputs",
                    {'fields':(("text", "html"),"summary"),
                     "classes": ("collapse",'extrapretty')}
                  ),

                 ("Run Time Settings",
                    {'fields':("json_text", "template"),
                     "classes": ("wide", 'extrapretty')},
                  ),
                 )


@admin.register(Project)
class ProjectAdmin(admin.ModelAdmin):

    formfield_overrides = {
        TextField: {'widget': Textarea(attrs={'rows':20,'cols': 100})},
    }

    search_fields = ('name', 'text', 'owner__first_name', 'owner__email')
    list_display = ("name", "date",  "deleted")
    list_filter = ( "name",  "deleted")

    fieldsets = (("Analysis Metadata",
                    {'fields': ("name","owner",("uid","sticky"),
                                 "deleted", "image"),
                     "classes": ('extrapretty')}
                  ),

                 ("Optional Text Inputs",
                    {'fields':("text","summary"),
                     "classes": ("collapse",'extrapretty')}
                  ),
                 )
