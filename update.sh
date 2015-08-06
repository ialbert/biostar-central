# The update script which is run when a githook payload is delivered. 
#
# This is a template script. Edit to your suit your deployment

# Kill the current serving process
waitressid=`ps aux | grep waitress |  awk -F ' ' '{print $2; exit;}'`
kill $waitressid

# Pull the changes
cd `pwd`
git pull

# Restart the server again
source conf/defaults.env
./biostar.sh waitress