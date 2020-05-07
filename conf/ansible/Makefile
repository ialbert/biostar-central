
# Absoulute path to the old biostar database on local machine
LOCAL_OLD_DB := data/old-biostar-database.sql.gz

# Github repository to use in the playbooks.
REPO := https://github.com/ialbert/biostar-central.git

# Target host must be set.
ifndef TARGET
$(error TARGET parameter must be set)
endif

setup:
	ansible-playbook -i hosts.ini -l ${TARGET} server-setup.yml

install:
	ansible-playbook -i hosts.ini -l ${TARGET} server-install.yml  --extra-vars "repo=${REPO}"

restart:
	ansible-playbook -i hosts.ini -l ${TARGET} server-restart.yml

transfer:
	echo LOCAL_OLD_DB=${LOCAL_OLD_DB}
	ansible-playbook -i ${HOST} server-migrate.yml --ask-become-pass --extra-vars "local_old_db=${LOCAL_OLD_DB} repo=${REPO}"  -v

deploy-code:
	ansible-playbook -i hosts.ini -l ${TARGET} server-deploy.yml  --extra-vars "repo=${REPO}"

deploy: deploy-code restart


# Temporary command to deploy forum at current state - with directory name being ~/biostar-engine/
# Being refactored out
forum_deploy:
	ansible-playbook -i hosts/test.biostars.org server-deploy.yml --ask-become-pass --extra-vars "runner=forum restart=True main_dir=/export/www/biostar-engine" -v
