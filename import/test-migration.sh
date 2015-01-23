#!/bin/bash
#
# Used to test federated content
#
rm -f biostar3/forum/migrations/0002_federatedcontent.py*
python manage.py makemigrations
./biostar.sh delete import init
python manage.py patch --federated_content import/federate/
./biostar.sh index
