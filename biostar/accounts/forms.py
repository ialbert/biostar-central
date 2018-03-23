from django import forms
from django.contrib.auth.models import User



class SignUpForm(forms.Form):

    password1 = forms.CharField(
        label="Password",
        strip=False,
        widget=forms.PasswordInput,
        max_length=254,
        min_length=2,
    )

    password2 = forms.CharField(
        label="Password confirmation",
        widget=forms.PasswordInput,
        strip=False,
        help_text="Enter the same password as before, for verification.",
    )

    email = forms.CharField(
        label="Email",
        strip=False,
        widget=forms.TextInput,
        max_length=254,
        min_length=2,
    )

    def clean_password2(self):

        password1 = self.cleaned_data.get("password1")
        password2 = self.cleaned_data.get("password2")
        if password1 and password2 and password1 != password2:
            raise forms.ValidationError("Passwords given do not match.")
        return password2

    def clean_email(self):

        data = self.cleaned_data['email']
        if User.objects.filter(email=data).exists():
            raise forms.ValidationError("This email is already being used.")
        return data


class LogoutForm(forms.Form):
    pass



class EditProfile(forms.Form):

    email = forms.CharField(label='Email', max_length=100)
    name = forms.CharField(label='First Name', max_length=100)

    def __init__(self, user,  *args, **kwargs):

        self.user = user

        super(EditProfile, self).__init__(*args, **kwargs)


    def save(self):

        self.user.email = self.cleaned_data['email']
        self.user.first_name = self.cleaned_data['name']
        self.user.save()

        return self.user

    def clean_email(self):

        data = self.cleaned_data['email']
        email = User.objects.exclude(pk=self.user.pk).filter(email=data)

        if email.exists():
            raise forms.ValidationError("This email is already being used.")

        return data


class LoginForm(forms.Form):

    email = forms.CharField(label='Email', max_length=100)
    password = forms.CharField(label='Password', max_length=100,
                               widget=forms.PasswordInput)
