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
        # give it a chance to load, not sure if this does what I think it does though
        for i in range(10):
            if part in text:
                break
            time.sleep(0.2)
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
    targets = "Latest,Latest questions,Follow multiple posts,Follow multiple tags,Follow multiple users".split(",")
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
    parts = "5 answers,Dondrup,uniprot,Biomart".split(",")
    contains(text, parts=parts)
    targets = "similar posts,permalink,revisions".split(",")
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
    fill(browser, "Other nearby island countries and territories include \
         Comoros and Mayotte to the south, Mauritius and Reunion to the far southeast,\
    and the Seychelles Islands about 1,500 km to the east", name='content')
    click(id="submit-button")
    click("zanzibar")
    
    # check the question appears on other pages
    click("Questions")
    click(title)
    
    # add an answer
    fill(browser, 'Take a boat then a plane then a train', name='content')
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
    feed_check,
    post_lookup,
    quick_search,
    update_user,
    create_content_1,
    voting_test,
]

def main(url):
    browser = webdriver.Firefox()
    browser.implicitly_wait(1)
    
    for func in tests:
        browser.get(url)
        func(browser)
    browser.close()

    
if __name__ == '__main__':
    import sys
    url="http://localhost:8080"
    
    if len(sys.argv)>1:
        url = sys.argv[1]
    
    main(url)