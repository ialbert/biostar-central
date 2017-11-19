from django.contrib import admin
from .models import Job, Analysis
from django.forms import TextInput, Textarea
from django.db.models import CharField, TextField



@admin.register(Job)
class JobAdmin(admin.ModelAdmin):

    formfield_overrides = {
        #CharField: {'widget': TextInput(attrs={'size': '20'})},
        TextField: {'widget': Textarea(attrs={'rows':20,'cols': 100})},
    }

    search_fields = ('name', 'owner__first_name', 'owner__email', 'state')
    list_display = ("name", "state","start_date", "security","date")
    list_filter = ("state", "security", "project__name")


    fieldsets = (("Analysis Metadata",
                    {'fields': ("name","owner",'project',("uid","sticky"),
                                ("state", "security")),
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

    search_fields = ('name', 'text', 'owner__first_name', 'owner__email' )
    list_display = ("name", "date", "auth")
    list_filter = ( "state", "auth", "project__name")


