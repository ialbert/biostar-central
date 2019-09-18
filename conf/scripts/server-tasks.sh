
 # * * * * *  command to execute
 # │ │ │ │ │
 # │ │ │ │ │
 # │ │ │ │ └───── day of week (0 - 6) (0 to 6 are Sunday to Saturday, or use names; 7 is Sunday, the same as 0)
 # │ │ │ └────────── month (1 - 12)
 # │ │ └─────────────── day of month (1 - 31)
 # │ └──────────────────── hour (0 - 23)
 # └───────────────────────── min (0 - 59)


# Send daily digests at 5am
@daily python digest.py --daily

# Send weekly digests
@weekly python digest.py --weekly

# Send monthly digest
@monthly  python digest.py --monthly


# Add posts to search index every 3 minutes
*/3 * * * *  python index.py