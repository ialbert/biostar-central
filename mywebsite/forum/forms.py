from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.utils.translation import gettext, gettext_lazy as helpers
from .models import User


class SignUpForm(UserCreationForm):

    #self = super()
    #super().__init__.self.instance.username =  super().__init__.self.cleaned_data.get('email')

    email = forms.EmailField(max_length=254, help_text='Required.')

    password1 = forms.CharField(max_length=254,
        min_length=2, 
        label=helpers("Password"),
        strip=False, 
        widget=forms.PasswordInput,
         )

    password2 = forms.CharField(
        label=helpers("Password confirmation"),
        widget=forms.PasswordInput,
        strip=False,
        help_text=helpers("Enter the same password as before."),
    )
    


    class Meta:
        model = User
        fields = ('email', 'password1', 'password2' )







