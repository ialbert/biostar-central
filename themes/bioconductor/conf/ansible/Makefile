

config:
	ansible-playbook -i ${HOST} server-config.yml -v

install:
	ansible-playbook -i ${HOST} server-install.yml --ask-become-pass --extra-vars -v

deploy:
	ansible-playbook -i ${HOST} server-deploy.yml --ask-become-pass --extra-vars "supervisor_program=${RUNNER} restart=True main_dir=${BASE_DIR}"  -v