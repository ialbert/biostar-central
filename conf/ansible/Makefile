
# Absoulute path to the old biostar database on local machine
LOCAL_OLD_DB := data/old-biostar-database.sql.gz

# Github repository to use in the playbooks.
REPO := https://github.com/ialbert/biostar-central.git

bioconductor:
	$(eval RUNNER:=forum)
	$(eval HOST := hosts/themes/supportupgrade.bioconductor.org)

config:
	ansible-playbook -i ${HOST} server-config.yml -v  --extra-vars "repo=${REPO} runner=${RUNNER}"

install:
	ansible-playbook -i ${HOST} server-install.yml -v --extra-vars "repo=${REPO} runner=${RUNNER}"

transfer:
	echo LOCAL_OLD_DB=${LOCAL_OLD_DB}
	ansible-playbook -i ${HOST} server-migrate.yml --ask-become-pass --extra-vars "local_old_db=${LOCAL_OLD_DB} repo=${REPO}"  -v

deploy:
	ansible-playbook -i ${HOST} server-deploy.yml --ask-become-pass --extra-vars "restart=True repo=${REPO} runner=${RUNNER}"  -v


# Temporary command to deploy forum at current state - with directory name being ~/biostar-engine/
# Being refactored out
forum_deploy:
	ansible-playbook -i hosts/test.biostars.org server-deploy.yml --ask-become-pass --extra-vars "runner=forum restart=True main_dir=/export/www/biostar-engine" -v
