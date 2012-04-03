from selenium import webdriver
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.common.keys import Keys
import time
from functools import partial
from django.conf import settings
from main.server import models

def click_func(browser, id, func=None):
    func = func or browser.find_element_by_link_text
    elem = func(id)
    assert elem, 'Element %s not found' % id
    elem.click()
    text = page(browser)
    return elem, text

def page(browser):
    return browser.page_source.encode("ascii", errors='replace')

def feed_check(browser):
    "Checks the feeds"
    click = partial(click_func, browser)
    elem, text  = click('About')
    elem, text  = click('Feeds page')
    targets = "Latest Questions,Follow multiple posts,Follow multiple tags,Follow multiple users".split(",")
    for link in targets:    
        elem, text = click(link) # Find the query box
        elem = browser.back()
        
def simple_navigation(browser):
    "Simple navigation through the site"
    click = partial(click_func, browser)
    targets = "Tags Users Badges About FAQ Recent Popular Questions Unanswered Planet Forum Tutorials Search!".split()
    targets.extend( [ 'My Tags', 'New Post!', 'Sign In' ] )
    for link in targets:    
        elem, text = click(link)

    drilldown = "Users,Istvan Albert,Bookmarks,Moderator,Tags,sequence".split(",")
    for link in drilldown: 
        elem, text = click(link)

def check(text, words):
    for word in words:
        assert word in text, "Unable to find %s in the text" % word

def fill(browser, name, text, submit=True):
    elem = browser.find_element_by_name(name)
    erase = Keys.BACK_SPACE * 100
    elem.send_keys(erase)
    submit = Keys.RETURN if submit else ''
    elem.send_keys(text + submit)
    time.sleep(0.5)
    text = page(browser)
    return text

def post_lookup(browser):
    click = partial(click_func, browser)
    link = "Gene ID conversion tool"
    elem, text = click(link)
    words = "5 answers,Dondrup,uniprot,Biomart".split(",")
    check(text, words=words)
    targets = "similar posts".split(",")
    for link in targets:
        click(link)
        browser.back()
    click("Renee")

def authenticate(uid, browser):
    click = partial(click_func, browser)
    link = "Sign In"
    
    settings.SELENIUM_TEST_LOGIN_TOKEN = 'murjkj4'
    settings.ACTIVE_USER = models.User.objects.get(id=uid)
    
    url = browser.current_url + "test/login/%s/%s/" % (uid, settings.SELENIUM_TEST_LOGIN_TOKEN)
    browser.get(url)
    text = page(browser)
    
    assert 'Test login ' in text

authenticate_user = partial(authenticate, 10)
authenticate_mod  = partial(authenticate, 2)

TITLE = "How to get to Zanzibar?"

def create_content(browser):
    click = partial(click_func, browser)
    click('New Post!')
    
    fill(browser, 'title', TITLE)
    fill(browser, 'tag_val', "zanzibar travel")
    fill(browser, 'content', "Other nearby island countries and territories include \
         Comoros and Mayotte to the south, Mauritius and Reunion to the far southeast,\
    and the Seychelles Islands about 1,500 km to the east")
    click("Submit Post")
    click("travel")
    
    # check the question appears on other pages
    click("Questions")
    click(tTITLE)
    
    # check the that search works
    fill(browser, 'q', "Mayotte")
    click(TITLE)

def add_answer(browser):
    click = partial(click_func, browser)
    
    # check the question appears on other pages
    click("Questions")
    click(TITLE)
    fill(browser, 'content', 'Take a boat then a plane then a train')
    
    
def update_user(browser):
    click = partial(click_func, browser)
    
    # modify the mytags settings
    user = settings.ACTIVE_USER
    click(user.profile.display_name)
    click("Edit info")
    fill(browser, 'scholar', "hordfUUAAAAJ", submit=False)
    fill(browser, 'my_tags', "travel")
    click("My Tags")
    time.sleep(2)
    
def detailed_navigation(browser):
    click = partial(click_func, browser)
    
def full_search(browser):
    "Searches via the main tab"
    click = partial(click_func, browser)
    click('Search!')
    
    elems = browser.find_elements_by_name('q')
    assert len(elems) == 2
    elem = elems[1]
    elem.send_keys("motif" + Keys.RETURN)
    time.sleep(0.2)
    text = page(browser)
    assert 'Finding common' in text
    
def quick_search(browser):    
    "Searches via the top box"
    
    text = fill(browser, 'q', "motif")
    assert 'Finding common' in text
    
    text = fill(browser, 'q', 'NO_SUCH_WORD_FOUND')
    assert 'found 0 results' in text
    
tests = [
    simple_navigation,
    detailed_navigation,
    post_lookup,
    quick_search,
    feed_check,
    authenticate_user,
    create_content,
    update_user,
    authenticate_mod,
    add_answer,
]

def main(url):
    browser = webdriver.Firefox()
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