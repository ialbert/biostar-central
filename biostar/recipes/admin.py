from django.contrib import admin
from django.db.models import TextField
from django.forms import Textarea

from .models import Job, Analysis, Project, Access, Data


@admin.register(Access)
class AccessAdmin(admin.ModelAdmin):
    list_select_related = (
        'user',
        'project',
    )

    readonly_fields = (
        'user',
    )

    pass


@admin.register(Data)
class DataAdmin(admin.ModelAdmin):
    pass


@admin.register(Job)
class JobAdmin(admin.ModelAdmin):
    formfield_overrides = {
        TextField: {'widget': Textarea(attrs={'rows': 20, 'cols': 100})},
    }

    search_fields = ('name', 'owner__first_name', 'owner__email', 'state', "project__name",
                     "project__owner__first_name", "project__owner__email")
    list_display = ("name", "state", "start_date", "security", "date")
    list_filter = ("state", "security", "project__name", "deleted")

    fieldsets = (("Job Metadata",
                  {'fields': ("name", "owner", 'project', ("uid"),
                              ("state", "security"), "image"),
                   "classes": ('extrapretty')}
                  ),

                 ("Optional Text Inputs",
                  {'fields': (("text", "html")),
                   "classes": ("collapse", 'extrapretty')}
                  ),

                 ("Run Time Settings",
                  {'fields': ("json_text", "template"),
                   "classes": ("wide", 'extrapretty')},
                  ),
                 )


@admin.register(Analysis)
class AnalysisAdmin(admin.ModelAdmin):
    formfield_overrides = {
        TextField: {'widget': Textarea(attrs={'rows': 20, 'cols': 100})},
    }

    search_fields = ('name', 'text', 'owner__first_name', 'owner__email', "project__name",
                     "project__owner__first_name", "project__owner__email")
    list_display = ("name", "date", "security")
    list_filter = ("security", "project__name", "deleted")

    fieldsets = (("Analysis Metadata",
                  {'fields': ("name", "owner", 'project', ("uid", "rank"),
                              ("deleted", "security"), "image"),
                   "classes": ('extrapretty')}
                  ),

                 ("Optional Text Inputs",
                  {'fields': (("text", "html")),
                   "classes": ("collapse", 'extrapretty')}
                  ),

                 ("Run Time Settings",
                  {'fields': ("json_text", "template"),
                   "classes": ("wide", 'extrapretty')},
                  ),
                 )


@admin.register(Project)
class ProjectAdmin(admin.ModelAdmin):
    formfield_overrides = {
        TextField: {'widget': Textarea(attrs={'rows': 20, 'cols': 100})},
    }

    search_fields = ('name', 'text', 'owner__first_name', 'owner__email')
    list_display = ("name", "date", "deleted")
    list_filter = ("name", "deleted")

    fieldsets = (("Analysis Metadata",
                  {'fields': ("name", "owner", ("uid", "rank"),
                              "deleted", "image", "privacy"),
                   "classes": ('extrapretty')}
                  ),

                 ("Optional Text Inputs",
                  {'fields': ("text",),
                   "classes": ("collapse", 'extrapretty')}
                  ),
                 )
