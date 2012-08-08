# Biostar version string

import sys, warnings
#warnings.simplefilter("default")
warnings.filterwarnings("ignore", category=DeprecationWarning) 

try:
    import django, south
except ImportError, exc:
    print "(!) unable to import a dependency"
    print "Try: unpack the depot.zip file in the libs folder"
    sys.exit()