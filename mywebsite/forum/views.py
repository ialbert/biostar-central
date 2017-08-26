from django.contrib.auth import login, authenticate
from django.shortcuts import render, redirect
from .forms import SignUpForm


def signup(request):

    if request.method == 'POST':
        form = SignUpForm(request.POST)
        # valid input at this point
        if form.is_valid():

            form.save()
            email = form.cleaned_data.get('email')

            raw_password = form.cleaned_data.get('password1')

            user = authenticate(email=email, password=raw_password)
            login(request, email)
            
            return redirect('/')
    else:
        
        form = SignUpForm()
    return render(request, 'forum/signup.html', {'form': form})



def login():
  return None



def index(request):
    return render(request, 'forum/home.html')





