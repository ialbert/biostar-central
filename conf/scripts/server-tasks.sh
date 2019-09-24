
 # * * * * *  command to execute
 # │ │ │ │ │
 # │ │ │ │ │
 # │ │ │ │ └───── day of week (0 - 6) (0 to 6 are Sunday to Saturday, or use names; 7 is Sunday, the same as 0)
 # │ │ │ └────────── month (1 - 12)
 # │ │ └─────────────── day of month (1 - 31)
 # │ └──────────────────── hour (0 - 23)
 # └───────────────────────── min (0 - 59)


# Send digests first minute of every day
@daily python ../manage.py digest.py --daily

# Send digests first minute of every week
@weekly python ../manage.py digest.py --weekly

# Send digests first minute of every month
@monthly  python ../manage.py digest.py --monthly

# Send watched tags emails first minute of every week
@weekly python ../manage.py digest.py --tags

# Add posts to search index every 3 minutes
*/3 * * * *  python ../manage.py index.py