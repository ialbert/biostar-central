from django import forms
from django.conf import settings
from django.core.exceptions import ValidationError
from main.server import const
import string, base64, json, hmac

P_TITLE, P_CONTENT, P_TAG = 'Post title', 'Post content', ''

# safe string transformation
import string
SAFE_TAG = set(string.ascii_letters + string.digits + "._- ")

def validate_integer(value):
    try:
        int(value)
    except (ValueError, TypeError), e:
        raise ValidationError('')

def valid_tag(text):
    "Validates form input for tags"

    if not text:
        raise ValidationError('Please enter at least one tag')

    if text == P_TAG:
        raise ValidationError('Please change the default tag')

    if len(text) > 300:
        raise ValidationError('Tag line is too long')

    words = text.split()

    if len(words) > 7:
        raise ValidationError('You have too many tags, please use at most seven tags')

def valid_title(text):
    "Validates form input for title"
    if text == P_TITLE:
        raise ValidationError('Please change the default title.')
    if len(text) < 5 :
        raise ValidationError('Your title appears to be shorter than the minimum of five characters.')
    if len(text) > 150 :
        raise ValidationError('Your title appears to be longer than the maximum of 150 characters.')

def valid_content(text):
    "Validates form input for content"
    # text size, min size, max size
    text = text.strip()
    ts, mi, mx = len(text), settings.MIN_POST_SIZE, settings.MAX_POST_SIZE
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

    title = forms.CharField(max_length=250, initial='', validators=[ valid_title ],
        widget=forms.TextInput(attrs={'class':'span8', 'placeholder': P_TITLE }))
    
    content = forms.CharField(max_length=10000, initial='', validators=[ valid_content ],
        widget=forms.Textarea(attrs={'cols':'80', 'rows':'15', 'id':'editor', 'placeholder': P_CONTENT}))

    tag_val = forms.CharField(max_length=250, initial='', validators=[ valid_tag ],
        widget=forms.TextInput(attrs={'class':'span4', 'placeholder': P_TAG}))

    context = forms.CharField(max_length=1000, required=False, initial='')

    type = forms.ChoiceField(choices=const.POST_TYPES[2:])

class ChildContent(forms.Form):
    """
    A form representing the body of simpler content answer/comment
    """
    content  = forms.CharField(max_length=10000, validators=[ valid_content ],
        widget=forms.Textarea(attrs={'cols':'80', 'rows':'15', 'id':'editor'}))

class RequestInfo(forms.Form):
    email = forms.CharField(max_length=200)

#
# Forms and validators for external login
#
def encode(data, key):
    text = json.dumps(data)
    text = base64.urlsafe_b64encode(text)
    digest = hmac.new(key, text).hexdigest()
    return text, digest

def decode(text, digest, key):
    if digest != hmac.new(key, text).hexdigest():
        raise Exception("message does not match the digest")
    text = base64.urlsafe_b64decode(text)
    data = json.loads(text)
    return data

class ExternalLogin(forms.Form):
    name = forms.CharField(max_length=50)
    data = forms.CharField(max_length=1500)
    digest = forms.CharField(max_length=64)
    action = forms.CharField(max_length=64, required=False)

    def clean(self):
        cleaned_data = super(ExternalLogin, self).clean()
        try:
            name = cleaned_data.get("name")
            data = str(cleaned_data.get("data"))
            digest = str(cleaned_data.get("digest"))
            if name not in settings.EXTERNAL_AUTHENICATION:
                raise Exception("unable to locate the authentication key by name")

            key, patt = settings.EXTERNAL_AUTHENICATION[name]
            data = decode(data, digest, key)
            username = data.get("username","").strip()
            if not username:
                raise Exception("username field not found in the data")
            cleaned_data['data'] = data
        except Exception, exc:
            raise forms.ValidationError("Invalid external login: %s" % exc)

        # Always return the full collection of cleaned data.
        return cleaned_data
