#!/bin/bash

# to be linked into .git/hooks/pre-push to fire off tests before push
# to bypass tests run
# git push --no-verify

CMD="python manage.py test --settings=run.sqlite"

# Check if we actually have commits to push
commits=`git log @{u}..`
if [ -z "$commits" ]; then
    exit 0
fi

$CMD
RESULT=$?
if [ $RESULT -ne 0 ]; then
	# Make it visible. Easy to miss.
	echo "***"
	echo "***"
    echo "*** Failed: $CMD"
    echo "***"
    echo "*** Push ERROR!"
    echo "***"
    echo "***"
    echo "***"
	# Try some beeps too.
	echo -ne '\a' '\a'
    exit 1
fi
exit 0
