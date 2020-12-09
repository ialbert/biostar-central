# Biostar API

This is the documentation for Biostar API. If you have additional questions, or believe you have encountered a bug, don't hesitate to post a question on Biostar.

## General
All API responses are JSON.

Some API responses are cached. Polling for changes should be done sparingly in any case, and polling at a rate faster than once a minute (for semantically identical requests) is considered abusive.

A number of methods in the Biostar API accept dates as parameters and return dates as properties, the format of these dates is documented above. As a general rule, full dates use ISO 8601 and timestamps are in unix epoch time.

## Methods
###Traffic

`GET /api/traffic/`

Number of post views over the last 60 min filtered by unique IPs.

#### Fields in response
- __date__: the current date, ISO 8601 format.
- __post_views_last_60_min__: number of post views over the last 60 min filtered by unique IPs.
- __timestamp__: the current date, unix epoch time format.

#### Example
/api/traffic/
```
{
    "date": "2014-05-29T14:59:55.788069",
    "post_views_last_60_min": 850,
    "timestamp": 1401375595
}
```
    
### User

`GET /api/user/{uid}/`

General info about a user.

#### Parameters
- __uid__: the identifier of the user.

#### Fields in response
- __date_joined__: the date the user joined the website, ISO 8601 format.
- __id__: the identifier of the user, a number.
- __joined_days_ago__: the date the user joined the website, as the number of days ago.
- __last_login__: the date of the last login of the user, ISO 8601 format.
- __name__: the name of the user.
- __vote_count__: the number of votes given by the user.

#### Example
/api/user/23/

```{
    "date_joined": "2010-01-18T21:43:55.253000+00:00",
    "id": 23,
    "joined_days_ago": 1614,
    "last_login": "2011-11-08T19:37:21.753000+00:00",
    "name": "Giovanni M Dall'Olio",
    "vote_count": 37
}
```
    
### Post

`GET /api/post/{id}/`

General info about a post.

#### Parameters
- __id__: the identifier of the post, a number.

#### Fields in response
- __answer_count__: number of answers.
- __author__: author name.
- __author_id__: author's identifier, a number.
- __book_count__: number of bookmarks.
- __comment_count__: number of comments.
- __creation_date__: creation date, ISO 8601 format.
- __has_accepted__: true if the question has an accepted answer, boolean.
- __id__: identifier of the post, a number.
- __lastedit_date__: date of last edit, ISO 8601 format.
- __lastedit_user_id__: user who last edited this post.
- __parent_id__: identifier of the parent post.
- __rank__: rank, a number.
- __reply_count__: number of replies.
- __root_id__: identifier of the root post.
- __status__: status message.
- __status_id__: status' identifier, a number.
- __subs_count__: number of subscribers following this post.
- __tag_val__: tags.
- __thread_score:__ thread's score.
- __title__: title.
- __type__: type of post.
- __type_id__: type's identifier for this post.
- __url__: url.
- __view_count__: number of views.
- __vote_count__: number of votes.
- __xhtml__: content.

Example
/api/post/25/
```
{
    "answer_count": 2,
    "author": "Gue Su",
    "author_id": 18,
    "book_count": 0,
    "comment_count": 0,
    "creation_date": "2009-12-01T20:57:35.300000+00:00",
    "has_accepted": false,
    "id": 25,
    "lastedit_date": "2009-12-01T20:57:35.300000+00:00",
    "lastedit_user_id": 18,
    "parent_id": 24,
    "rank": 0.0,
    "reply_count": 0,
    "root_id": 24,
    "status": "Open",
    "status_id": 1,
    "subs_count": 0,
    "tag_val": "",
    "thread_score": 0,
    "title": "A: How To Set Shrimp Parameters For Best Sensitivity With 35Bp Colorspace Data?",
    "type": "Answer",
    "type_id": 1,
    "url": "http://localhost:8080/p/24/#25",
    "view_count": 0,
    "vote_count": 2,
    "xhtml": "
I just read the SHRiMP manual again, but I think that their explanation about -M option may not be enough to answer your question. I usually use the \"seed\" mode by using -s, -n, and -w and the option -M is a new feature of the version 1.3.1, which I have never tried before.

\n\n
I recommend for you to use the \"seed\" mode--the default would be good, but please adjust the -s option if you want more sensitivity. Always fast speed compensates sensitivity and the -M option seems to exist for this purpose.

\n\n
Hope my message to be helpful for your project.

\n"
}
```

### Vote

`GET /api/vote/{id}/`

General info about a vote.

#### Parameters
- __id__: the identifier of the vote, a number.

#### Fields in response
- __author__: author name.
- __author_id__: author's identifier, a number.
- __date__: date of the vote, ISO 8601 format.
- __id__: identifier of the vote, a number.
- __post_id__: identifier of the voted post.
- __type__: type of vote.
- __type_id__: type's identifier for this vote.

Example
/api/vote/21/
```
{
    "author": "Zhaorong",
    "author_id": 14,
    "date": "2014-04-29T15:02:17.740000+00:00",
    "id": 21,
    "post_id": 26,
    "type": "Upvote",
    "type_id": 0
}
```
    
### Statistics on the Nth day

`GET /api/stats/day/{day}/`

Statistics as of the Nth day after day-0 (the day of the first ever post).

#### Parameters
- __day__: number of days after day-0, a number.

#### Fields in response
- __answers__ : total number of answers as of the given day.
- __comments__: total number of comments as of the given day.
- __date__: date, ISO 8601 format.
- __new_posts__: number of new posts in the given day.
- __new_users__: number of new users in the given day.
- __new_votes__: number of new votes in the given day.
- __questions__: total number of questions as of the given day.
- __timestamp__: date, unix epoch time format.
- __toplevel__: total number of toplevel post as of the given day.
- __users__: total number of users as of the given day.
- __votes__: total number of votes as of the given day.

#### Example
/api/stats/day/5/
```
{
    "answers": 6,
    "comments": 0,
    "date": "2009-10-05T00:00:00",
    "new_posts": [
        10,
        11,
        12
    ],
    "new_users": [
        10,
        11
    ],
    "new_votes": [],
    "questions": 6,
    "timestamp": 1254700800,
    "toplevel": 6,
    "users": 10,
    "votes": 0
}
```
    
### Statistics on a date

`GET /api/stats/date/{year}/{month}/{day}/`

Statistics as of the given date.

#### Parameters
- __year__: a number, 4 digits.
- __month__: a number, 2 digits.
- __day__: a number, 2 digits.

#### Fields in response
- __answers__: total number of answers as of the given date.
- __comments__: total number of comments as of the given date.
- __date__: date, ISO 8601 format.
- __new_posts__: number of new posts in the given date.
- __new_users__: number of new users in the given date.
- __new_votes__: number of new votes in the given date.
- __questions__: total number of questions as of the given date.
- __timestamp__: date, unix epoch time format.
- __toplevel__: total number of toplevel post as of the given date.
- __users__: total number of users as of the given date.
- __votes__: total number of votes as of the given date.

Example
/api/stats/date/2009/10/06/

```
{
    "answers": 9,
    "comments": 0,
    "date": "2009-10-06T00:00:00",
    "new_posts": [
        13,
        14,
        15,
        16
    ],
    "new_users": [
        12,
        13
    ],
    "new_votes": [],
    "questions": 7,
    "timestamp": 1254787200,
    "toplevel": 7,
    "users": 12,
    "votes": 0
}
  ```
  
### Tags List

`POST /api/tags/list/`

Return a list of tags with corresponding counts of posts. Can also pass down a time range.

#### Parameters
- __data__: a file listing the tags, with 
- __months__: 6


Given 

    curl -X POST -F "tags=@/Users/natay/Desktop/apps/biostar-central/tags.txt" http://localhost:8000/api/tags/list/?trange=year
    

Returns 

    {
        "tag1": {
            "answer_count": 0,
            "comment_count": 0,
            "total": 21
        },
        "tag2": {
            "answer_count": 0,
            "comment_count": 0,
            "total": 20
        }
    }