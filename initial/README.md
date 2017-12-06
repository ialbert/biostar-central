# Server install

sudo apt-get update && sudo apt-get upgrade -y

apt-get install -y fail2ban

sudo apt-get install -y curl unzip build-essential ncurses-dev
sudo apt-get install -y byacc zlib1g-dev python-dev git cmake
sudo apt-get install -y screen parallel
sudo apt-get install -y default-jdk ant
sudo apt-get install -y postgresql nginx

# useradd -m -s /bin/bash www

# mkdir -p /export


source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("DESeq")
biocLite("DESeq2")
biocLite("edgeR")
install.packages("gplots")