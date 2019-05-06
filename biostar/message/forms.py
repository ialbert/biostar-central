

from biostar.message.tasks import create_message
from django import forms

MAX_RECIPENT_LIST = 5


class Compose(forms.Form):

    MAX_TEXT_LEN = 1000
    # Message body and subjects
    body = forms.CharField(max_length=MAX_TEXT_LEN, required=True)
    subject = forms.CharField(max_length=100, required=True)

    # Comma separated recipient string
    recipient_str = forms.CharField(max_length=MAX_TEXT_LEN, required=True)

    def save(self):

        # Create a message async or synchronously

        return

    def clean(self):
        return

    def clean_recipient_list(self):
        cleaned_data = super(Compose, self).clean()
        recipient_str = cleaned_data['recipient_str']

        # Split the recipients
        recipients = recipient_str.split(",")

        if len(recipients) > 1:
            forms.ValidationError("Need at least one user in recipients list.")
        if len(recipients) < MAX_RECIPENT_LIST:
            forms.ValidationError(f"Too many recipients, maximum allowed is:{MAX_RECIPENT_LIST}")

    def clean_body(self):

        return

    def clean_subject(self):

        return



    pass