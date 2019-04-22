## Deploying the Biostar Engine

The site is built with Django therefore the official Django documentation applies to maintaining and deploying the site:

* <https://docs.djangoproject.com/>

## Running jobs

A recipe submitted for execution is called a job.

When the job is run the recipe parameters are applied onto recipe template to produce the script that gets executed. This transformation takes place right before executing the job.

Jobs can be executed as commands. See the `job` command for details:

    python manage.py job --help

The command has number of parameters that facilitate job management and recipe development.
For example:

    python manage.py job --list

will list all the jobs in the system. Other flags that allow users to investigate and override the behaviors.

    python manage.py job --id 4 --show_script

will print the script for job 4 that is to be executed to the command line. Other flags such as `-use_template` and `-use_json` allows users to override the data or template loaded into the job.
This can be useful when developing new recipes.

Another handy command:

    python manage.py job --next

will execute the next queued job. The job runner may be run periodically with cron.

## Automatic job spooling

The Biostar Engine supports `uwsgi`. When deployed through
`uwsgi` jobs are queued and run automatically through the `uwsgi` spooler. See the `uwsgi` documentation  for details on how to control and customize that process.

* <https://uwsgi-docs.readthedocs.io/en/latest/>

[uwsgi]: <https://uwsgi-docs.readthedocs.io/en/latest/

## Security consideration

**Note**: The site is designed to execute scripts on a remote server. In addition the site
allows users with **moderator** rights may change the content of these scripts.

It is **extremely important** to monitor, restrict and guard access to all
accounts with moderator privileges!
