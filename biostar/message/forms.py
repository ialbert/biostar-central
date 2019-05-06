
from django.forms import forms


class Compose(forms.Form):

    # Message content and subjects
    body = ""
    subject = ""

    # Comma separated recipient list
    recipient_list = ""


    def clean(self):
        return

    pass