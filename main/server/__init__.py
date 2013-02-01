import sys

# this is the BioStar release number
VERSION = '1.2.15'

try:
    import docutils
except ImportError, exc:
    print '(!) unable to import the docutils module'
    print '(!) some functionality may be disabled'
    #sys.exit(-1)