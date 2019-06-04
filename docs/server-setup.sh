# Must be run as root.

set -uex

# The host name of the server.
HOST_NAME=test.biostars.org

# Set the hostname.
hostnamectl set-hostname $HOST_NAME

# Set the timezone.
dpkg-reconfigure tzdata

# Add the default user.
adduser www

# Add the www user to the sudo group.
adduser www sudo

# Update the server.
apt-get update && apt-get upgrade -y

# Upgrade the distributions.
sudo apt-get dist-upgrade

# Add and enable fail2ban.
apt-get install ufw fail2ban -y

# Enable the firewalls
ufw allow ssh
ufw allow http
ufw allow https
ufw enable

# Set up a packages that will be needed.
apt-get install nginx postgresql software-properties-common curl -y
apt-get install supervisor build-essential cmake zlib1g-dev -y

# Install the certbot packages.
add-apt-repository universe
add-apt-repository ppa:certbot/certbot
apt-get update
apt-get install certbot python-certbot-nginx -y

