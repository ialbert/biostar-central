from django import forms

P_TITLE, P_CONTENT, P_TAG = 'Post title', 'Post content', 'tag1'
class PostForm(forms.Form):
    """
    A form representing a new question
    """
    title   = forms.CharField(max_length=250,  initial=P_TITLE)
    content = forms.CharField(max_length=5000, initial=P_CONTENT)
    tags    = forms.CharField(max_length=250,  initial=P_TAG)
    
    def clean(self):
        "Custom validator for the question"
        if not super(PostForm, self).is_valid():
            raise forms.ValidationError("Invalid form")
        
        if self.cleaned_data['tags'] == P_TAG:
            raise forms.ValidationError("Please create a different tag")
    
        if self.cleaned_data['content'] == P_CONTENT:
            raise forms.ValidationError("Please change the content")
        
        if self.cleaned_data['title'] == P_TITLE:
            raise forms.ValidationError("Please create a different title")

        return self.cleaned_data
    
class ContentForm(forms.Form):
    """
    A form representing the body of simpler content answer/comment
    """
    content = forms.CharField(max_length=5000)

