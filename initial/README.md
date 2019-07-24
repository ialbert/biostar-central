# Various other server installation

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
biocLite("DESeq", ask=FALSE)
biocLite("DESeq2",ask=FALSE)
biocLite("edgeR", ask=FALSE)
install.packages("gplots", quiet=TRUE)


conda create -y --name py2 python=2
source activate py2
pip install fabric
    
curl https://img.shields.io/badge/biostar-engine-red.svg?style=flat-square > biostar/engine/static/images/badge.svg
curl https://img.shields.io/badge/biostar-recipes-green.svg?style=flat-square > biostar/engine/static/images/recipe-badge.svg

