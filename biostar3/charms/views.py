from __future__ import absolute_import, division, print_function, unicode_literals
from django.shortcuts import render, redirect
import json, time
from django.http import HttpResponse
from decorator import decorator

@decorator
def json_response(func, request):
    """
    Returns data as a http response in JSON format
    """
    return HttpResponse(json.dumps(func(request)))

def charms(request):
    """
    This view creates nodes. Is not called directly from the web only through
    other functions that prefill parameters.
    """
    template_name = "charm_view.html"
    context = dict()
    return render(request, template_name, context)

@json_response
def charms_rpc(request):

    #time.sleep(3)

    command = request.POST.get("command")

    result = "{}".format(command)

    status = "error"
    msg = "Job request submitted. Data will appear later on this page."
    result = "1"

    data = dict(status=status, msg=msg, command=command, result=result)

    print (data)

    return data
