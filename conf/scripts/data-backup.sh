PATT=`date  +"%Y-%m-%d"`

SITE=http://localhost:8000

TOKEN='92d9b420-d912-40'
LIST=$SITE/api/list/
PROJ=$SITE/api/project/

mkdir -p data

curl $LIST?token=$TOKEN > recipes.json

cat recipes.json | jq -r 'keys[]' | parallel -j 1 "curl ${PROJ}{}/?token=${TOKEN} > data/{}.json"

tar czvf recipes-json-${PATT}.tar.gz recipes.json data/*