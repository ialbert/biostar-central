PATT=`date  +"%Y-%m-%d"`

SITE=http://localhost:8000

LIST=$SITE/api/list/
PROJ=$SITE/api/project/

mkdir -p data

curl $LIST > recipes.json

cat recipes.json | jq -r 'keys[]' | parallel -j 1 "curl ${PROJ}{}/ > data/{}.json"

tar czvf recipes-json-${PATT}.tar.gz recipes.json data/*