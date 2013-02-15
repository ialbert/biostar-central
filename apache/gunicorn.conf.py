bind = 'unix:/tmp/gunicorn.sock'

#bind = 'localhost:8001'

workers = 2
debug = True
max_requests = 1000
loglevel="debug"

print "*** loading gunicorn configuration"
