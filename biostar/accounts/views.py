import uuid, logging

from django.contrib import messages
from ratelimit.decorators import ratelimit
from django.shortcuts import render, redirect
from django.urls import reverse
from django.contrib.auth.models import User
from django.contrib.auth import views as auth_views
from django.utils.safestring import mark_safe

from .forms import SignUpForm, LoginForm, LogoutForm, EditProfile
from .models import Profile
from django.contrib import auth

logger = logging.getLogger('engine')


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def edit_profile(request):

    if request.user.is_anonymous:
        messages.error(request, "Must be logged in to edit profile")
        return redirect("/")

    id = request.user.id
    user = User.objects.filter(id=id).first()

    if request.method == "POST":
        form = EditProfile(data=request.POST, user=user)
        if form.is_valid():
            form.save()
            return redirect("profile")

        messages.error(request, mark_safe(form.errors))

    initial = dict(email=user.email, first_name=user.first_name)
    form = EditProfile(initial=initial)
    context = dict(user=user, form=form)
    return render(request, 'accounts/edit_profile.html', context)


def profile(request):

    if request.user.is_anonymous:
        messages.error(request, "Must be logged in to edit profile")
        return redirect("/")

    id = request.user.id
    user = User.objects.filter(id=id).first()
    context = dict(user=user)

    return render(request, 'accounts/profile.html', context)


@ratelimit(key='ip', rate='10/m', block=True, method=ratelimit.UNSAFE)
def user_signup(request):

    messages.info(request, "Signups are not yet enabled")
    return redirect("/")

    # steps = breadcrumb_builder([HOME_ICON, SIGNUP_ICON])
    #
    # if request.method == 'POST':
    #
    #     form = SignUpForm(request.POST)
    #     form.add_error(None, "Sign up is disabled")
    #     if form.is_valid():
    #         email = form.cleaned_data.get('email')
    #         password = form.cleaned_data.get('password1')
    #         name = email.split("@")[0]
    #
    #         user = User.objects.create(username=get_uuid(), email=email,
    #                                    first_name=name)
    #         user.set_password(password)
    #         user.save()
    #
    #         auth.login(request, user)
    #         logger.info(f"Signed up and logged in user.id={user.id}, user.email={user.email}")
    #         messages.info(request, "Signup successful!")
    #         return redirect(reverse('login'))
    # else:
    #     form = SignUpForm()
    # context = dict(form=form, steps=steps)
    # return render(request, 'accounts/signup.html', context=context)


def user_logout(request):

    if request.method == "POST":

        form = LogoutForm(request.POST)

        if form.is_valid():
            auth.logout(request)
            messages.info(request, "You have been logged out")
            return redirect("/")

    form = LogoutForm()

    context = dict(form=form)

    return render(request, "accounts/logout.html", context=context)


def user_login(request):

    if request.method == "POST":
        form = LoginForm(data=request.POST)

        if form.is_valid():

            email = form.cleaned_data['email']
            password = form.cleaned_data['password']

            # Due to an early bug emails may not be unique. Last subscription wins.
            user = User.objects.filter(email__iexact=email).order_by('-id').first()

            if not user:
                messages.error(request, "This email does not exist.")
                return redirect(reverse("login"))

            user = auth.authenticate(username=user.username, password=password)

            if user and user.profile.state in [Profile.BANNED, Profile.SUSPENDED ]:
                msg = f"Login not allowed. Account is <b> {user.profile.get_state_display()}</b>"
                messages.error(request, mark_safe(msg))
                return redirect("/")

            if not user:
                messages.error(request, "Invalid Password")
                return redirect(reverse("login"))

            elif (user and not user.is_active):
                messages.error(request, "This user may not log in.")
                return redirect(reverse("login"))

            elif user and user.is_active:
                auth.login(request, user)
                messages.success(request, "Login successful!")
                return redirect("/")


        messages.error(request, mark_safe(form.errors))

    form = LoginForm()
    context = dict(form=form)
    return render(request, "accounts/login.html", context=context)



def password_reset(request):
    context = dict()

    return auth_views.password_reset(request, extra_context=context,
                                     template_name="accounts/password_reset_form.html",
                                     subject_template_name="accounts/password_reset_subject.txt",
                                     email_template_name="accounts/password_reset_email.html"
                                     )

def password_reset_done(request):
    context = dict()

    return auth_views.password_reset_done(request, extra_context=context,
                                          template_name="accounts/password_reset_done.html")

def pass_reset_confirm(request, uidb64, token):
    context = dict()

    return auth_views.password_reset_confirm(request, extra_context=context,
                                             template_name="accounts/password_reset_confirm.html",
                                            uidb64=uidb64, token=token)

def password_reset_complete(request):
    context = dict()

    return auth_views.password_reset_complete(request, extra_context=context,
                                              template_name="accounts/password_reset_complete.html",)
