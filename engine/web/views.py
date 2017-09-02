from django.shortcuts import render


def index(request):
    return render(request, 'index.html')

def blog(request):
    return render(request, 'blog/blog.html')

def data(request):
    return render(request, 'blog/data.html')

def forum(request):
    return render(request, 'blog/forum.html')
