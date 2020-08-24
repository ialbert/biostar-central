
# Absoulute path to the old biostar database on local machine
LOCAL_OLD_DB := data/old-biostar-database.sql.gz

# Github repository to use in the playbooks.
REPO := https://github.com/ialbert/biostar-central.git

USER := root

# Target host must be set.
ifndef TARGET
$(error TARGET parameter must be set)
endif

setup:
	ansible-playbook -i hosts.ini -l ${TARGET} server-setup.yml

install:
	ansible-playbook -i hosts.ini -l ${TARGET} server-install.yml  --extra-vars "repo=${REPO} username=${USER}"

restart:
	ansible-playbook -i hosts.ini -l ${TARGET} server-restart.yml --extra-vars "repo=${REPO} username=${USER}"

transfer:
	echo LOCAL_OLD_DB=${LOCAL_OLD_DB}
	ansible-playbook -i ${HOST} server-migrate.yml --ask-become-pass --extra-vars "local_old_db=${LOCAL_OLD_DB} repo=${REPO} username=${USER}"  -v

deploy-code:
	ansible-playbook -i hosts.ini -l ${TARGET} server-deploy.yml  --extra-vars "repo=${REPO} username=${USER}"

deploy: deploy-code restart


# Temporary command to deploy forum on test server - with directory name being ~/biostar-engine/
# Being refactored out
forum_deploy:
	ansible-playbook -i hosts.ini -l test server-deploy.yml --extra-vars "main_dir=/export/www/biostar-engine username=${USER} repo=${REPO}" -v
	ansible-playbook -i hosts.ini -l test server-restart.yml --extra-vars "main_dir=/export/www/biostar-engine username=${USER} repo=${REPO}"
