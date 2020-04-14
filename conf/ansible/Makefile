
RUNNER := engine

recipes:
	$(eval HOST := hosts/www.bioinformatics.recipes)

config:
	ansible-playbook -i ${HOST} server-config.yml -v

install:
	ansible-playbook -i ${HOST} server-install.yml -v

deploy:
	ansible-playbook -i ${HOST} server-deploy.yml --ask-become-pass --extra-vars "restart=True"  -v

# Temporary command to deploy forum at current state - with directory name being ~/biostar-engine/
# Being refactored out

forum_deploy:
	ansible-playbook -i ${HOST} server-deploy.yml --ask-become-pass --extra-vars "supervisor_program=forum restart=True main_dir=/export/www/biostar-engine" -v
