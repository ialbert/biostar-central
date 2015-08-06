from django.shortcuts import render, redirect, render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.views.decorators.csrf import csrf_exempt
from django.utils.crypto import constant_time_compare
import json, hmac, base64, hashlib
from netaddr import IPAddress, IPNetwork
from django.conf import settings
import os, json

def verifyIP(ip):
	if IPAddress(ip) in IPNetwork('192.30.252.0/22'):
		return True
	else:
		return False

def verifyUser(payload):
	jsonbody = json.loads(payload)
	pusher = jsonbody['pusher']['name']
	handle = settings.GITHUB_HANDLE
	if(pusher == handle):
		return True
	else:
		return False


@csrf_exempt
def webhookupdate(request):
    if request.META.get('CONTENT_TYPE') == 'application/json':
    	payload = request.body
        if verifyIP(request.META.get('REMOTE_ADDR')) and verifyUser(payload):
            cmd = settings.HOME_DIR+"/update.sh"
	    os.system('%s'%(cmd))
	    return HttpResponse("All done!")
        else:
            raise Http404
    else:
        raise Http404
