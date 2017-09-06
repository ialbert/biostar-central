from django import forms
from django.utils.translation import gettext, gettext_lazy as helpers
from django.contrib.auth.models import User
from django.contrib.auth.forms import AuthenticationForm


class SignUpForm(forms.ModelForm):
    
    password1 = forms.CharField(
        label=helpers("Password"),
        strip=False,
        widget=forms.PasswordInput,
        max_length=254,
        min_length=2,
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


    def clean_password2(self):

        password1 = self.cleaned_data.get("password1")
        password2 = self.cleaned_data.get("password2")
        if password1 and password2 and password1 != password2:
            raise forms.ValidationError(
                helpers("Passwords given do not match."))
        return password2

    def clean_email(self):

        data = self.cleaned_data['email']
        if User.objects.filter(email=data).exists():
            raise forms.ValidationError("This email is already being used.")
        return data

    def save(self, commit=True):

        user = super().save(commit=False)
        
        if commit:
            user.save()
            return
        return user 


    def cleaned_data(self, *args):
        return self.cleaned_data


class CustomLoginForm(AuthenticationForm):
    def confirm_login_allowed(self, user):
        print( user.is_validated)
        if not user.is_active or not user.is_validated:
            raise forms.ValidationError('There was a problem with your login.', 
                        code='invalid_login')









