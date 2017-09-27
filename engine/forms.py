from django import forms
from django.utils.translation import gettext_lazy as helpers
from django.contrib.auth.models import User
from .models import Project, Data, Analysis, Result
from django.forms.extras.widgets import SelectDateWidget

from pagedown.widgets import PagedownWidget

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

    def cleaned_data(self, *args):
        return self.cleaned_data


class LoginForm(forms.Form):

    email = forms.CharField(label='Email', max_length=100)
    password = forms.CharField(label='Password', max_length=100, 
                              widget=forms.PasswordInput)


class ProjectForm(forms.ModelForm):

    text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))

    class Meta:

        model = Project
        fields = ['title', 'text']


class DataForm(forms.ModelForm):

    text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))
    class Meta:
        model = Data
        fields = ['title', 'text']


    def cleaned_data(self, *args):
        return self.cleaned_data


class AnalysisForm(forms.ModelForm):

    text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))
    class Meta:
        model = Analysis
        fields = ['title', 'text']


class RunForm(forms.Form):
    # reads from json file here?

    QUEUED, RUNNING, FINISHED, ERROR = 1,2,3,4

    CHOICES = [(QUEUED, "Queued"), (RUNNING, "Running"),
               (FINISHED, "Finished"), (ERROR, "Error")]

    COMMANDS = [(1, "setup_dir"),
                (2, "check_input"),
                (3, "create_multiqc")]

    state = forms.IntegerField(label="State", widget=forms.Select(choices=CHOICES) )
    commands = forms.IntegerField(label="Commands", widget=forms.Select(choices=COMMANDS))

    def __init__(self, *args, **kwargs):

        analysis_id = kwargs.pop("id")
        super().__init__(*args, **kwargs)

        analysis = Analysis.objects.filter(id=analysis_id).first()
        data = analysis.project.data_set.all()

        CHOICES = []

        for i,d in enumerate(data):
            CHOICES.append((i, d))

        self.fields["data"] = forms.CharField(label="Pick Data",
                                              widget=forms.Select(choices=CHOICES) )

class ResultForm(forms.ModelForm):

    text = forms.CharField(widget=PagedownWidget(template="widgets/pagedownwidget.html"))

    class Meta:
        model = Result
        fields = ['title', 'text', 'commands', 'state', 'directory']




