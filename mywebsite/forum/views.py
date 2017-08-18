from django.shortcuts import render


def index(request):
    return render(request, 'forum/header.html')

