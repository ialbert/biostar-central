# Scripts available to the recipes app


## api_backup.sh

Run api back up. Edit `TOKEN` and `SITE` 

One place where requests are made in parallel:

    # Download each project into its own file
    cat recipes.json | jq -r 'keys[]' | parallel -j 1 "curl ${PROJ}{}/?token=${TOKEN} > data/{}.json"
    
## api_upload.sh


Run api upload on zipped file found in `DUMP`. This file is returned from `api_backup.sh`.

Edit `TOKEN` and `SITE` .

    # 
    ls data/ |  parallel -j 1 "curl -X POST -F 'data=@{}' ${UPLOAD} "