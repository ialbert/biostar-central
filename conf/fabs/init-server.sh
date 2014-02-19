#!/bin/sh
#
# Initializing an Ubuntu based linux distribution
#
IP="0.0.0.0"
HOSTNAME="hostname"

# Set the hostname
echo  $HOSTNAME > /etc/hostname
hostname -F /etc/hostname

# Initialize the hosts file
echo "127.0.0.1        localhost.localdomain    localhost" > /etc/hosts
echo "$IP      $HOSTNAME        $HOSTNAME" >> /etc/hosts

# reconfigure timezone
dpkg-reconfigure tzdata

# Update linux distro to latest
apt-get update -y
apt-get upgrade -y --show-upgraded

# Install postgresql
apt-get install -y postgresql nginx fail2ban redis-server ufw

# Start installing required packages
apt-get install -y build-essential ncurses-dev byacc zlib1g-dev python-dev git
apt-get install -y python-setuptools

# Install pip
easy_install pip

# Install the virtual environments
pip install virtualenv virtualenvwrapper

# Set up and start ufw
ufw allow ssh
ufw allow http
ufw enable

# Start services
service nginx start