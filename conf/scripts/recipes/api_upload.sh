
SITE=http://localhost:8000


# Secure API token, match user profile token
TOKEN='change-me'

DUMP=recipes-json-latest.tar.gz
UPLOAD=$SITE/api/upload/

# Specify what file to send
#UID=demo

tar -zxvf ${DUMP}

ls data/ |  parallel -j 1 "curl -X POST -F 'data=@{}' ${UPLOAD}/?token=${TOKEN} "
