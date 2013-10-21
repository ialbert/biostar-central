from django.conf import settings

from celery import task
from django.core.mail import send_mail
from django.contrib.sites.models import Site

@task
def demo_add(x, y):
    return x + y

@task
def send_test_email():
    subject = "test email"
    sender  = "istvan.albert@gmail.com"
    recipient = "iua1@psu.edu"
    body = "this is a test email"
    send_mail(subject, body, sender, [recipient], fail_silently=False)

if __name__ == '__main__':
    send_test_email()
