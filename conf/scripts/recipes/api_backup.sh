PATT=`date  +"%Y-%m-%d"`

SITE=http://localhost:8000

# Secure API token, match user profile token
TOKEN='change-me'

LIST=$SITE/api/list/
PROJ=$SITE/api/project/

# Make directory to collect files
mkdir -p data

# Fetch list of available projects
curl $LIST?token=$TOKEN > recipes.json

# Download each project into its own file
cat recipes.json | jq -r 'keys[]' | parallel -j 1 "curl ${PROJ}{}/?token=${TOKEN} > data/{}.json"

# Zip into one archive.
tar czvf recipes-json-${PATT}.tar.gz recipes.json data/*