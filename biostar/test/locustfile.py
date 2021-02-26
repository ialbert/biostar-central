import time
from locust import HttpUser, task, between

class QuickstartUser(HttpUser):
    wait_time = between(1, 2.5)

    @task
    def home(self):
        self.client.get("/")

    #@task
    def visit_page(self):
        self.client.get("/p/2/")

    #@task
    def visit_popular(self):
        self.client.get("/p/7126/")

