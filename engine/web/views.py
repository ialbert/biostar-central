from random import choice
from string import ascii_lowercase, digits
import uuid
from django.contrib import auth
from django.contrib.auth import login, authenticate
from django.contrib.auth import views as auth_views
from django.shortcuts import render, redirect
from .forms import SignUpForm
from django.contrib.auth import get_user_model
from ratelimit.decorators import ratelimit

class LoginForm(forms.Form):
    email = forms.CharField(label='Email', max_length=100)
    password = forms.CharField(label='Password', max_length=100, widget=forms.PasswordInput)

def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


@ratelimit(key='ip', rate='10/m', block=True, method=ratelimit.UNSAFE)
def login(request):
    if request.method == "POST":
        form = forms.LoginForm(request.POST)

        if form.is_valid():
            email = form.cleaned_data['email']
            password = form.cleaned_data['password']

            # Due to an early bug emails may not be unique. Last subscription wins.
            user = User.objects.filter(email__iexact=email).order_by('-id').first()

            if not user:
                form.add_error(None, "This email does not exist.")
                context = dict(form=form)
                return render(request, "registration/user_login.html", context=context)

            user = auth.authenticate(username=user.username, password=password)

            if not user:
                form.add_error(None, "Invalid password.")
            elif user and not user.is_active:
                form.add_error(None, "This user may not log in.")
            elif user and user.is_active:
                auth.login(request, user)
                logger.info("logged in user.id={}, user.email={}".format(user.id, user.email))
                return redirect("/")
            else:
                # This should not happen normally.
                form.add_error(None, "Invalid form processing.")
    else:
        initial = dict(nexturl=request.GET.get('next', '/'))
        form = forms.LoginForm(initial=initial)

    context = dict(form=form)
    return render(request, "registration/user_login.html", context=context)



def already_used(username):

    usermodel = get_user_model()
    try:
        usermodel.objects.get(username=username)
        return True
    except usermodel.DoesNotExist:
        return False


def hash_username(length=16):

    chars=ascii_lowercase+digits
    username = ''
    for token in range(length):
        username += choice(chars)

    if already_used(username):
        hash_username()
    else:
        return username 


def signup(request):

    if request.method == 'POST':
        form = SignUpForm(request.POST)

        if form.is_valid():

            email = form.cleaned_data.get('email')
            raw_password = form.cleaned_data.get('password1')

            user = User.objects.create(username=hash_username(), email=)
            user = form.save(commit=False)
            user.username = hash_username()
            user.set_password(raw_password)
            user.email = email
            user.save()

            login(request, user)
            
            return redirect('/login')
    else:
        
        form = SignUpForm()
    return render(request, 'signup.html', {'form': form})



def index(request):
    return render(request, 'index.html')

   


