
RUNNER := engine

# Absoulute path to the old biostar database on local machine
LOCAL_OLD_DB := data/old-biostar-database.sql.gz


recipes:
	$(eval HOST := hosts/www.bioinformatics.recipes)

config:
	ansible-playbook -i ${HOST} server-config.yml -v

install:
	ansible-playbook -i ${HOST} server-install.yml -v

transfer:
	echo LOCAL_OLD_DB=${LOCAL_OLD_DB}
	ansible-playbook -i ${HOST} server-migrate.yml --ask-become-pass --extra-vars "local_old_db=${LOCAL_OLD_DB}"  -v

deploy:
	ansible-playbook -i ${HOST} server-deploy.yml --ask-become-pass --extra-vars "restart=True"  -v


# Temporary command to deploy forum at current state - with directory name being ~/biostar-engine/
# Being refactored out

forum_deploy:
	ansible-playbook -i hosts/test.biostars.org server-deploy.yml --ask-become-pass --extra-vars "supervisor_program=forum restart=True main_dir=/export/www/biostar-engine" -v
