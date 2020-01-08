#!/bin/bash

#
# tested local crontab with 
# */10 * * * * /bin/bash --login -c 'cd /home/lori/SupportSiteTagsScript/generateTags && ./supportSiteTags.sh > /home/lori/SupportSiteTagsScript/generateTags/GenerateTags.log 2>&1'
#

# Add additional specialized tags
echo -e "Bioconductor\nnew package\nsingle package builder\nsubmission\nrelease\ndevel\ninstall\ndeprecation\ndeadline\nconference\nworkshop\nbuild report\nreports\nwebsite" > tags.txt

# software packages
curl https://bioconductor.org/packages/devel/bioc/VIEWS | grep "Package:" | cut -d ' ' -f 2  >> tags.txt

# experiment data packages
curl https://bioconductor.org/packages/devel/data/experiment/VIEWS | grep "Package:" | cut -d ' ' -f 2  >> tags.txt 

# annotation packages
curl https://bioconductor.org/packages/devel/data/annotation/VIEWS | grep "Package:" | cut -d ' ' -f 2  >> tags.txt

# workflow packages
curl https://bioconductor.org/packages/devel/workflows/VIEWS | grep "Package:" | cut -d ' ' -f 2  >> tags.txt

# biocViews vocab terms
curl https://raw.githubusercontent.com/Bioconductor/biocViews/master/inst/dot/biocViewsVocab.dot | grep "\->" | sed '/^\//d' | sed -r 's/[A-Za-z]*\s->\s([A-Za-z]*);/\1/'  >> tags.txt

# remove misc carriage return
sed -i 's/\r//g' tags.txt

# sort and remove duplicates
sort -u -o tags.txt tags.txt


# check for curl failure
checkFile=`grep "curl:" GenerateTags.log | wc -l`
fileDiff=`diff ../tags.txt tags.txt | grep '<\|>' | wc -l`

if [ "$checkFile" -gt "0" ]
then
    echo "Possible curl failure. Not updating tags.txt"
elif [ "$fileDiff" -eq "0" ]
then
    echo "Files do not differ. Not updating tags.txt"
else
    echo "Files differ. Attempting to update tags.txt"
    cp tags.txt ../
fi
