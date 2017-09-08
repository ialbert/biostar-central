from abc import abstractmethod, ABCMeta
from django.db import models


class A(models.Model):

    a=3
    b=4



class B(A):
    print(a)



d = B()

print(d)
