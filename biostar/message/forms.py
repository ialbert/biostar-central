

from biostar.message.tasks import send_message
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
        subject = self.cleaned_data.get("subject", "subject")
        body = self.cleaned_data.get("body", "Body")
        rec_str = self.cleaned_data("recipient_str")
        rec_list = rec_str.split(",")
        # Use the IN query on a known lwngth
        rec_users = User.objects.filter(usernae__in=rec_str.split(","))
        # Create a message
        send_message(subject, body, rec_list, sender)
        return

    def clean(self):
        return

    def clean_recipient_str(self):
        cleaned_data = super(Compose, self).clean()
        recipient_str = cleaned_data['recipient_str']

        # Split the recipients
        recipients = recipient_str.split(",")

        if len(recipients) > 1:
            forms.ValidationError("Need at least one user in recipients list.")
        if len(recipients) < MAX_RECIPENT_LIST:
            forms.ValidationError(f"Too many recipients, maximum allowed is:{MAX_RECIPENT_LIST}")
        return recipient_str

    def clean_body(self):

        return

    def clean_subject(self):

        return



    pass