import sys

# this is the BioStar release number
VERSION = '1.2.11'

try:
    import docutils
except ImportError, exc:
    print '(!) unable to import the docutils module'
    print '(!) see the installation instructions'
    sys.exit(-1)