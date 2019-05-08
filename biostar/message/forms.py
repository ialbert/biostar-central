from pagedown.widgets import PagedownWidget
from biostar.accounts.models import User
from biostar.message.tasks import send_message
from django import forms

MAX_RECIPIENT_LIST = 5


class Compose(forms.Form):

    MAX_TEXT_LEN = 5000
    # Message body and subjects
    body = forms.CharField(widget=PagedownWidget(template="widgets/pagedown.html"),
                           min_length=2, max_length=MAX_TEXT_LEN,
                           help_text="Message body.")

    subject = forms.CharField(max_length=100, required=True, help_text="Message subject.")

    # Comma separated recipient string
    recipients = forms.CharField(max_length=100, required=True,
                                 help_text="Comma separated list of user handles: user1, user2, etc...")

    def save(self, sender):
        subject = self.cleaned_data.get("subject", "subject")
        body = self.cleaned_data.get("body", "Body")
        rec_str = self.cleaned_data.get("recipients", "")
        rec_list = rec_str.split(",")
        # Use the __in query on a known length list
        rec_users = User.objects.filter(username__in=rec_list).exclude(username=sender)
        # Create a message
        send_message(subject=subject, body=body, rec_list=rec_users, sender=sender)
        return

    def clean(self):
        return

    def clean_recipients(self):
        cleaned_data = super(Compose, self).clean()
        recipients = cleaned_data['recipients']

        # Split the recipients
        rec_list = recipients.split(",")

        if len(rec_list) < 1:
            forms.ValidationError("Need at least one user in recipients list.")
        if len(rec_list) < MAX_RECIPIENT_LIST:
            forms.ValidationError(f"Too many recipients, maximum allowed is:{MAX_RECIPIENT_LIST}")
        return recipients