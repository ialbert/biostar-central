
# Frequently Asked Questions

<a name="top"></a>

This is the FAQ for the *Bioconductor* Support Site. Please visit
[Bioconductor.org](https://bioconductor.org/) for additional project and package information. A list of
current *Bioconductor* packages can be found [here](https://bioconductor.org/packages/release/BiocViews.html)

  * [General](#general)
    * [What is the difference between the support site and the bioc-devel mailing list?](#g1)
    * [What questions should I post?](#g2)
    * [How do I contact you?](#g3)
    * [What if I find a bug with the support site?](#g4)
    * [Where can I suggest new features to the support site?](#g5)
    * [What about offensive content?](#g6)
  * [Posts](#posts)
    * [How do I write good posts?](#p1)
    * [How do I format a post?](#p2)
    * [How do I put images into my posts?](#p3)
    * [What post tags are automatically suggested?](#p4)
    * [How can I keep track of topics that interest me?](#p5)
    * [What are the different post types?](#p6)
    * [How do I post a Job, Tutorial, or News item?](#p7)
    * [How do I filter posts?](#p8)
  * [Account](#account)
    * [How do I merge multiple previous acccounts?](#a1)
    * [User reputation](#a2)
    * [How do I edit my profile?](#a3)
    * [How do I update my profile picture?](#a4)
    * [What is a handle and how do I update it?](#a5)
    * [What are the different notification options?](#a6)
    * [What is the different between 'My Tags' and 'Watched Tags'?](#a7)
    * [Does selecting a digest option change my notification option?](#a8)
    * [How are awards generated?](#a9)
    * [What awards are available to earn?](#a10)
    * [Becoming a moderator](#a11)
  * [Banner](#banner)
    * [What shows up in 'Messages'?](#b1)
    * [What shows up in 'Votes'?](#b2)
    * [What shows up in 'My Posts'](#b3)
    * [What shows up in 'My Tags'](#b4)
    * [What shows up in 'Follow'?](#b5)
    * [What shows up in 'Bookmarks'? / How do I bookmark a post?](#b6)
    * [What does selecting a digest option do? Does it effect my notification option?](#b7)
  * [Other](#other)

<a name="general"></a>

# General

<a name="g1"></a>

#### What is the difference between the support site and the bioc-devel mailing list?

Do not post to both locations! The support site is for questions related to
usage of *Bioconductor* packages and troubleshooting; This site is appropriate
for most questions. For guidance on developing a *Bioconductor* package,
including current and future implementation of functions and features, please
use the bioc-devel mailing list

<a name="g2"></a>

#### What questions should I post?

The support site is for questions related to usage of *Bioconductor* packages.
*Bioconductor* package maintainers are required upon submission to subscribe
to the support site and are expected to answer questions regarding their
packages at this location. General questions about bioinformatics,
computational genomics, biological data analysis and software unrelated to
Bioconductor can be asked on [biostars](https://www.biostars.org/) or
[seqanswers](http://seqanswers.com/). General R questions can be asked on the
[R-help mailing list](https://stat.ethz.ch/mailman/listinfo/r-help) or
[StackOverflow](https://stackoverflow.com/questions/tagged/r).

<a name="g3"></a>

#### How do I contact you?

Those of us who maintain this site read the posts that are made here. So users
who wish to contact us are encouraged to post to this site.

<a name="g4"></a>

#### What if I find a bug with the support site?

No site is perfect and you may come across things that you feel are not
working as they should be. This site is based on a fork of Biostars so we have
the power to correct many problems as they are discovered. If you find a
problem with this website, we ask that you please post such questions to our
[github repository](https://github.com/Bioconductor/support.bioconductor.org/issues).

<a name="g5"></a>

#### Where can I suggest new features to the support site?

New feature requests can be made by opening an issue on the [github
repository](https://github.com/Bioconductor/support.bioconductor.org/issues).

<a name="g6"></a>

#### What about offensive content?

Users posting content that does not belong to the site will be notified and
required to edit their content. Users may post job postings to the 'Jobs'
section as long as the topic aligns with the main focus of this site. Users
posting obvious spam will be immediately suspended or banned. Moderators will
try to remove offensive content in a timely fashion but if you find something
offensvie please notify webmaster@bioconductor.org. Post not relating to
*Bioconductor* content will generally be removed by moderators.

<div style="text-align: right"> [ <a href="#top">Back to top</a> ] </div>

<a name="posts"></a>

# Posts

<a name="p1"></a>

#### How do I write good posts?

One of the biggest things that new users struggle with is how to ask questions
in a way that allows people to answer them. So if you are unfamiliar please
have a look at the [posting
guide](http://www.bioconductor.org/help/support/posting-guide/) for some great
advice and the [tutorial](https://support.bioconductor.org/p/117436/) for help
on using the markdown editor. If a *Bioconductor* package reports an error that
you cannot debug and are requesting assistance, Please check and include the
following:

  * That `BiocManager::valid()` reports your packages are of the correct version.
  * Record the output of `traceback()` immediately after the error occurs.
  * Include the output of `sessionInfo()`.
  * Provide minimal code necessary to try and reproduce the issue. StackOverflow has suggestions for making a [reproducible example](https://stackoverflow.com/help/minimal-reproducible-example).

<a name="p2"></a>

#### How do I format a post?

The support site usese a markdown editor. There is a useful tutorial post with
additional helpful links found [here](https://support.bioconductor.org/p/117436/).

<a name="p3"></a>

#### How do I put images into my posts?

You can also put images and plots into your posts. All you need to do is
create (or login to) a free account with http://imgur.com and upload any
images or plots that you want to share. Then you can enter the imgur generated
'Direct Link' url into the 'URL' location that is available when you click on
the 'image' button in the text editor. The plot will be imbedded right into
your post.

<a name="p4"></a>

#### What post tags are automatically suggested?

When creating a post, a user is required to enter at least one 'Post Tag'.
There will be automatically suggested tags. These tags are generated from the
current list of *Bioconductor* packages as well as biocViews.

<a name="p5"></a>

#### How can I keep track of topics that interest me?

This site offers two ways to track posts when they use certain tags. You can
add the tags that you want tracked into your user profile. 'My Tags' will put
posts that match your tagging criteria into a special tab that you can check
whenever you log in. And 'Watched Tags' should actually email you whenever
someone posts using a tag that you are interested in.

<a name="p6"></a>

#### What are the different post types?

When ASK A QUESTION is selected there are four options. Most often a user
should keep the "Ask a question" as the selected 'Post Type'. The four options
are:

  * Ask a question: For general questions and help
  * Post a Job Opening: For recruiting new hires by posting open job opportunities for a related field appropriate to the *Bioconductor* community
  * Share a Tutorial: Share a tutorial for a bioinformatic analysis or *Bioconductor* package
  * Announce News: Announce breaking news about *Bioconductor* or package to the community at large

<a name="p7"></a>

#### How do I post a Job, Tutorial, or News item?

Follow the ASK A QUESTION link and choose the appropriate field in the 'Post
Type' drop down box.

<a name="p8"></a>

#### How do I filter posts?

There are a number of different ways to filter posts. A user can edit their
'My Tags' in their profile to allow specific posts in the 'My Tags' section of
the banner. To filter the main posting list, a user could search for a
specific term in the search field above the post listings. A user could also
filter using the 'Limit' and 'Sort' drop down boxes above the post listings.

<div style="text-align: right"> [ <a href="#top">Back to top</a> ] </div>

<a name="account"></a>

# Account

<a name="a1"></a>

#### How do I merge multiple previous acccounts?

So perhaps you have been using the Bioconductor mailing list for a long time
and during that time you have been using multiple different aliases. And now
you would like to consolidate things so that you can get credit for actually
being only one person. What should you do? Well 1st of all you should get set
up with the account that you normally use. If you have not logged in before
this may mean going straight to password recovery.

After setting up the primary account, that contact the site admins to resolve the issue. 

<a name="a2"></a>

#### User reputation

The number next to a user's name is the sum of upvotes and accepted answers
that user has collected. Upvotes boost your reputation on the site.

<a name="a3"></a>

#### How do I edit my profile?

Navigate to your user account by either clicking on your name or go to the
list of users, search for your user name, and click on it. From this page, you
should be able to select 'Edit profile' under your gravatar/profile picture.

<a name="a4"></a>

#### How do I update my profile picture?

<a name="a5"></a>

#### What is a handle and how do I update it?

A handle is a convenient way for other users to tag you in a post so you
receive a email notification regarding the posting. If a users handle was
'user1', anyone could tag the user in a post with '@user1'. This triggers an
email notification for the post. A user can change their handle by editing
their profile. See previous section on editing profile.

<a name="a6"></a>

#### What are the different notification options?

When editing your profile there are different notification options to choose
from. The following describe the options:

  * Default:
  * Email:
  * Local Messages:
  * No Messages:

<a name="a7"></a>

#### What is the different between 'My Tags' and 'Watched Tags'?

'My Tags' controls which posts can be seen in the 'My Tags' tab. Any post with
the tags listed will show up in the tab. 'Watched Tags' controls email
postings. Any post with these tags will trigger an email notification about
the post. It is expected that package maintainers have their package selected
for at minimum 'Watched Tags'.

<a name="a8"></a>

#### Does selecting a digest option change my notification option?

<a name="a9"></a>

#### How are awards generated?

Awards are triggered automatically by performing actions or reaching
milestones required for the award.

<a name="a10"></a>

#### What awards are available to earn?

<a name="a11"></a>

#### Becoming a moderator

Moderators are selected by *Bioconductor* core team members. The
responsibility of moderating posts is given to long-term trusted users.

<div style="text-align: right"> [ <a href="#top">Back to top</a> ] </div>

<a name="banner"></a>

# Banner

<a name="b1"></a>

#### What shows up in 'Messages'?

The Messages tab shows awards earned and comments on posts.

<a name="b2"></a>

#### What shows up in 'Votes'?

The Votes tab shows posts that have received an upvote, accepted answer, or if
someone has bookmarked your post.

<a name="b3"></a>

#### What shows up in 'My Posts'

Any post that you have created, answered or commented appear in this tab.

<a name="b4"></a>

#### What shows up in 'My Tags'

Any post with the tags identified in 'My Tags' of a user's profile are listed
under this tab. These can be edited by updating/editing your user profile.

<a name="b5"></a>

#### What shows up in 'Follow'?

<a name="b6"></a>

#### What shows up in 'Bookmarks'? / How do I bookmark a post?

After selecting a post, there is a bookmark icon underneath the upvote count
for the main post, as well as any answer, comment, or reply.

<a name="b7"></a>

#### What does selecting a digest option do? Does it effect my notification option?

<div style="text-align: right"> [ <a href="#top">Back to top</a> ] </div>

<a name="other"></a>


####  What HTML tags allowed in the post content

Font styling:
   -  `b`, `i`,`img`, `strong`, `strike`, `em`, `underline`
   - `sub`, `sup`, `img`
   - `h1` , `h2`, `h3`, `h4`
   
Content layouts:
   - `div`, `span`, `br`,`p`, `hr`
   
Code render
   - `code`, `pre` 
   
Tables
   -`table`, `thead`, `tr`, `th`, `td`, `tbody`

Allowed styles:

   - `color`
   - `font-weight`
   - `background-color`
   - `width height`
   
# Other

Can't find what your looking for? Consider opening an issue on the
support.bioconductor.org [github issues](https://github.com/Bioconductor/support.bioconductor.org/issues).

<div style="text-align: right"> [ <a href="#top">Back to top</a> ] </div>
