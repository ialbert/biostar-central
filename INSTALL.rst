To run the biostar Q&A site you will need to have Python and Django 1.3 installed.

To check for python type::

	python --version

It needs to be Python 2.6 or higher.

To check if you have the necessary Django version, run:
  ./run.sh init

To install Django 1.3 without affecting your hosts system-wide settings,
you can use python-virtualenv so you will not need to use root.

1. This is the only time root access is needed. Install virtualenv::

	sudo apt-get install python-virtualenv
	# or
    sudo easy_install virtualenv

2. Create a python environment folder where Django will be installed::

	virtualenv env

3. Activate the environment::

	source env/bin/activate 

4. Install Django 1.3 without becoming root::
	
	pip install django

5. Now you should be able to run biostar-central
   cd ../biostar-central
  ./run.sh init
