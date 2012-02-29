from django import forms
from django.core.exceptions import ValidationError
from main.server import const
import string

P_TITLE, P_CONTENT, P_TAG = 'Post title', 'Post content', 'tag1'

def validate_integer(value):
    try:
        int(value)
    except (ValueError, TypeError), e:
        raise ValidationError('')

def valid_tag(text):
    "Validates form input for tags"
    
    if not(text):
        raise ValidationError('Please enter at least one tag')

    if text == P_TAG:
        raise ValidationError('Please change the default tag')
    
    words = text.split()
    
    if len(words) > 5:
        raise ValidationError('You have too many tags, please use at most five tags')
    
    for word in words:
        if len(word) < 3:
            raise ValidationError("Tag '%s' is too short, use at least 3 characters" % word)
        if len(word) > 16:
            raise ValidationError("Tag '%s' is too long, use no more than 16 characters" % word)

def valid_title(text):
    "Validates form input for title"
    if text == P_TITLE:
        raise ValidationError('Please change the default title.')
    if len(text) < 5 :
        raise ValidationError('Your title appears to be shorter than the minimum of five characters.')
    if len(text) > 100 :
        raise ValidationError('Your title appears to be longer than the maximum of 100 characters.')

def valid_content(text):
    "Validates form input for content"
    # text size, min size, max size
    text = text.strip()
    ts, mi, mx = len(text), 15, 10000
    if not(text):
        raise ValidationError('Content appears to be whitespace')
    if text == P_CONTENT:
        raise ValidationError('Please change the default content')
    if ts < mi :
        raise ValidationError('Your content is only %d charactes long. The minimum is %d.' %(ts, mi))
    if ts > mx :
        raise ValidationError('Your content  is too long %d characters. The maximum is %d .' % (ts, mx))
  
class TopLevelContent(forms.Form):
    """
    A form representing a new question
    """
    error_css_class = 'error'
    required_css_class = 'required'

    title = forms.CharField(max_length=250,  initial=P_TITLE, validators=[ valid_title ],
        widget=forms.TextInput(attrs={'class':'span8', 'onfocus':"remove(this, '%s')" % P_TITLE }))
    
    content = forms.CharField(max_length=10000, initial=P_CONTENT, validators=[ valid_content ], 
        widget=forms.Textarea(attrs={'cols':'80', 'rows':'15', 'id':'editor', 'onfocus':"remove(this, '%s')" % P_CONTENT}))

    tag_val = forms.CharField(max_length=250,  initial=P_TAG, validators=[ valid_tag ], 
        widget=forms.TextInput(attrs={'class':'span4', 'onfocus':"remove(this, '%s')" % P_TAG}))
    
    # the first two post types are not creatable here
    type = forms.ChoiceField(choices=const.POST_TYPES[2:])

class ChildContent(forms.Form):
    """
    A form representing the body of simpler content answer/comment
    """
    content  = forms.CharField(max_length=5000, validators=[ valid_content ],
        widget=forms.Textarea(attrs={'cols':'80', 'rows':'15', 'id':'editor'}))


