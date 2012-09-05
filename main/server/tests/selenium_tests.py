from selenium import webdriver
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.common.keys import Keys
import time
from functools import partial
from django.conf import settings
from main.server import models, const

def page(browser):
    return browser.page_source.encode("ascii", errors='replace')

def contains(text, parts):
    for part in parts:
        assert part in text, "Unable to find %s in the text" % part
        
def get(browser, link=None, name=None, id=None):
    if link:
        print '*** link:%s' % link
        return browser.find_element_by_partial_link_text(link)
    elif name:
        print '*** name:%s' % name
        return browser.find_element_by_name(name)
    elif id:
        print '*** id:%s' % id
        return browser.find_element_by_id(id)
    else:
        raise Exception("no identifiers were set")
        
def click_func(browser, link=None, name=None, id=None):
    elem = get(browser, name=name, id=id, link=link)
    assert elem, 'Element name=%s, id=%s not found' % (name, id)
    elem.click()
    text = page(browser)
    return elem, text

def feed_check(browser):
    "Checks the feeds"
    click = partial(click_func, browser)
    
    elem, text  = click('About')
    elem, text  = click('Feeds page')
    targets = "New Posts,New Questions,Follow multiple posts,Follow multiple tags,Follow multiple users".split(",")
    for link in targets:    
        elem, text = click(link)
        elem = browser.back()
        
def simple_navigation(browser):
    "Simple navigation through the site"
    click = partial(click_func, browser)
    
    targets = "Tags Users Badges About Recent Planet Search Posts News Questions Unanswered Tutorials Tools Videos Jobs RSS".split()
    targets.extend( [ 'Posts', 'Show All', 'New Post!',  'Sign In' ] )
    for link in targets:    
        elem, text = click(link)

    drilldown = "Users,Istvan Albert,Bookmarks,Moderator,Tags,sequence,Badges,Teacher".split(",")
    for link in drilldown: 
        elem, text = click(link)


def fill(browser, text, name=None, id=None, submit=False):
    elem  = get(browser, name=name, id=id)
    erase = Keys.BACK_SPACE * 100
    submit = Keys.RETURN if submit else ''
    elem.send_keys(erase + text + submit)
    time.sleep(0.5)
    text = page(browser)
    return text
    
def post_lookup(browser):
    click = partial(click_func, browser)
    
    click("next>")
    target = "Gene ID conversion tool"
    elem, text = click(target)
    time.sleep(2)
    parts = "Dondrup,uniprot,Biomart,Agilent".split(",")
    contains(text, parts=parts)
    targets = "similar posts,link,revisions".split(",")
    for link in targets:
        click(link)
        browser.back()
        
    # the username that created this post
    click("Renee")
    
    click("Planet")
    #blog = models.Post.objects.filter(type=const.POST_BLOG)[0]
    #click(blog.title)
    #browser.back()
    
def login(browser, uid):
    "Logs in with a user"
    print '*** logging in user %s' % uid
    user = models.User.objects.get(id=uid)
    url = browser.current_url + "test/login/%s/%s/" % (uid, settings.SELENIUM_TEST_LOGIN_TOKEN)
    browser.get(url)
    text = page(browser)
    contains(text, [ "login complete" ])
    return user

def create_content_1(browser):    
    # logging in as a general user

    root = browser.current_url
    user = login(browser=browser, uid=10)
    
    click = partial(click_func, browser)

    click('New Post!')

    title = "How to get to Zanzibar?"
    fill(browser, title, name='title')
    fill(browser, "zanzibar", name='tag_val' )
    content = """
    
Other nearby island countries and territories include 
Comoros and Mayotte to the south, Mauritius and Reunion to the far southeast,
and the Seychelles Islands about 1,500 km to the east

Orphan links should be autolinked: http://www.biostars.org, same with ftp://www.biostars.org

Links within codeblocks should be kept verbatim:

    http://www.biostars.org
    
We have enabled a number of new features on Biostar. Using these is optional and
are only meant to facilitate certain use cases.

Shortcuts
=========

Shortcuts are words that start with a backslash followed by one or more
comma separated parameters.

For example one if a user would like to share
code via [Gist][gist] they may write:

    \gist 2059

This shortcut will be replaced inside the post by \gist 2059

Shortcut properties:

- more than one comma separated parameter value may be listed
- shortcuts will not take effect if they are in a codebox (as in the example).
- unrecognized shortcuts will pass into the main text unchanged.
- the preview will not show the results of the shortcuts, users need to submit
the content to have them take effect.

[gist]: https://gist.github.com/

Shorcuts for Users and Posts
=============================

Shortcuts to link to users and posts:

    \user 23
    \post 34
    
Results in \user 23 and \post 34

For each of the examples you may list multiple values separated by commas:

    \user 23,30
    \post 34,79
    
Results in \user 23,30 and \post 34,79

Shorcuts for Tags
=================

    \tag blast,pipeline
    
Results in the link \tag blast,pipeline

Embedding Gist
==============

[Gist][gist] is a simple way to share snippets and pastes with others.
All gists are git repositories, so they are automatically versioned, forkable and usable as a git repository 

    \gist 2059
    
Result: \gist 2059

Embed Searchbox
===============

To create a filled in search box write:

    \search this
    
Result: \search this

Embed Youtube
=============

Use the Youtube ID as below (you can extract this ID from the url of the video):

    \youtube 1ZyoI-4ObSA

Results in \youtube 1ZyoI-4ObSA
    """
    fill(browser, content, name='content')
    click(id="submit-button")
    click("zanzibar")
    
    # check the question appears on other pages
    click("Questions")
    click(title)
    answer="""
Take a boat then a plane then a train
But let's also test orpan linking here: http://www.biostars.org but not insided
code blocks:

    http://www.biostars.org

First link should be a real html link, second should stay as is.
"""
    # add an answer
    fill(browser, answer, name='content')
    click(id="submit-button")
    
def update_user(browser):
    click = partial(click_func, browser)
    root = browser.current_url
    
    user = login(browser=browser, uid=10)
    click(user.profile.display_name)
    click("Edit info")
    fill(browser, "mapping", name="my_tags")
    fill(browser, "Cool Bot", name="display_name")
    click(id='submit-button')
    click("Posts")
    click("My Tags")    
    title = "Gene ID conversion tool"
    click(title)
 
    click("RSS")
    click("Activity Feed")
    browser.back()
    click("My Tags Feed")

def full_search(browser):
    "Searches via the main tab"
    click = click_link(browser)
    click('Search')
    
    elems = browser.find_elements_by_name('q')
    assert len(elems) == 2
    elem = elems[1]
    elem.send_keys("motif" + Keys.RETURN)
    time.sleep(0.2)
    text = page(browser)
    assert 'Finding common' in text
    
def quick_search(browser):    
    "Searches via the top box"
    
    click = partial(click_func, browser)
    
    text = fill(browser, "motif", name='q', submit=True)
    contains(text, parts=[ "Finding common motifs",  ])
    
    text = fill(browser, 'NO_SUCH_WORD_FOUND', name='q', submit=True)
    contains(text, parts=[ 'found 0 results' ])
  
def voting_test(browser):
    title = "Gene ID conversion tool"
    
    user = login(browser=browser, uid=10)
    
    click = partial(click_func, browser)
    title = "How to organize a pipeline of small scripts together?"
    click(title)
    
    post = models.Post.objects.get(title=title)
    author = post.author
    votes  = list(models.Vote.objects.filter(author=user, post=post, type=const.VOTE_UP))
    
    get_score = lambda x: models.User.objects.get(id=x).profile.score
    # get the scores
    
    score1 = get_score(author.id)
    first  = browser.find_elements_by_class_name("vote-up")[0]
    first.click()
    time.sleep(1)
    
    score2 = get_score(author.id)
    diff = score2 - score1
    
    # existing vote, removes vote
    if votes:
        assert diff== -1, "Voting error 1"
    else:
        assert diff == +1, "Voting error 2"

    # second round of voting
    votes  = list(models.Vote.objects.filter(author=user, post=post, type=const.VOTE_UP))
    first.click()
    time.sleep(1)
    score3 = get_score(author.id)
    diff = score3 - score2
    
    if votes:
        assert diff== -1, "Voting error 1"
    else:
        assert diff == +1, "Voting error 2"
    
    # it is already bookmarked it will remove the bookmark
    bookmarked = models.Vote.objects.filter(author=user, post=post, type=const.VOTE_BOOKMARK).count()
    
    # check bookmarking
    bookmark  = browser.find_elements_by_class_name("vote-bookmark")[0]
    bookmark.click()
    time.sleep(2)
    count = models.Vote.objects.filter(author=user, post=post, type=const.VOTE_BOOKMARK).count()
    if bookmarked:
        assert count == 0, 'Bookmark has not been removed'
    else:
        assert count == 1, 'Item was not bookmarked'
        click(user.profile.display_name)
        click('Bookmarks')
        click(title)
        
       
    
tests = [
    simple_navigation,
    post_lookup,
    quick_search,
    update_user,
    create_content_1,
    voting_test,
    feed_check,
]

#tests = [ create_content_1 ]

def main(url):
    browser = webdriver.Firefox()
    browser.implicitly_wait(1)
    
    for func in tests:
        browser.get(url)
        # give a chance to the page to load up
        # not sure why this is needed
        time.sleep(2)
        func(browser)
    browser.close()

    
if __name__ == '__main__':
    import sys
    url="http://localhost:8080"
    
    if len(sys.argv)>1:
        url = sys.argv[1]
    
    main(url)