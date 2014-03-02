#!/bin/sh
# Initializes the filesystem
set -ue

echo "IP=$IP"
echo "HOSTNAME=$HOSTNAME"

BIOSTAR_HOME=~/sites/$HOSTNAME

mkdir -p ~/.ssh
mkdir -p ~/backup
mkdir -p $BIOSTAR_HOME

# User needs to check the above for validity.
read -p "Press [Enter] key to start server prep..."

cd $BIOSTAR_HOME



