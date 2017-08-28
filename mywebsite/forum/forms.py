from django import forms
from django.contrib.auth import password_validation
from django.utils.translation import gettext, gettext_lazy as helpers
from .models import User


class SignUpForm(forms.ModelForm):
 
    error_messages = {
        'password_mismatch': helpers("The two password fields didn't match."),
    }
    password1 = forms.CharField(
        label=helpers("Password"),
        strip=False,
        widget=forms.PasswordInput,
    )
    password2 = forms.CharField(
        label=helpers("Password confirmation"),
        widget=forms.PasswordInput,
        strip=False,
        help_text=helpers("Enter the same password as before, for verification."),
    )

    class Meta:

        model = User
        fields = ("email",)

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        if self._meta.model.USERNAME_FIELD in self.fields:
            self.fields[self._meta.model.USERNAME_FIELD].widget.attrs.update({'autofocus': True})

    def clean_password2(self):

        password1 = self.cleaned_data.get("password1")
        password2 = self.cleaned_data.get("password2")
        if password1 and password2 and password1 != password2:
            raise forms.ValidationError(
                helpers("The two password fields didn't match."))
        return password2

    def save(self, commit=True):

        user = super().save(commit=False)
        user.set_password(self.cleaned_data["password1"])
        if commit:
            user.save()
        return user

    def cleaned_data(self, *args):
        return self.cleaned_data





