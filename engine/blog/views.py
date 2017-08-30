from django.shortcuts import render


def index(request):
    return render(request, 'home.html')

def blog(request):
    return render(request, 'blog.html')

def data(request):
    return render(request, 'data.html')


def forum(request):
    return render(request, 'forum.html')
