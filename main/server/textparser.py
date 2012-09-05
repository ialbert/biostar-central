import sys
import ply
import ply.lex as lex

class MainLexer:
    # List of token names.   This is always required
    states = (
        ('pre', 'exclusive'),
        ('post','exclusive'),
    )

    tokens = (
       'NEWLINE',
       'SPACE',
       'LINK',
       'REFERENCE',
       'ORPHAN',
       'CODE_OPEN', 'CODE_CLOSE',
       'TAG_OPEN', 'TAG_CLOSE',
       'COMMAND',
       'TEXT',
    )

    def t_pre_NEWLINE(self, t):
        r'\n|\r\n'
        self.lexer.code_block = False
        self.lexer.start_line = True
        return t

    def t_post_NEWLINE(self, t):
        r'\n|\r\n'
        self.lexer.start_line = True
        return t

    def t_pre_SPACE(self, t):
        r'(\s+)'
        if len(t.value) == 4 and self.lexer.start_line:
            self.lexer.code_block = True
        return t

    def t_post_SPACE(self, t):
        r'(\s+)'
        return t

    def t_pre_REFERENCE(self, t):
        r'\[\s*\S+\s*\]\:\s*\S+'
        self.lexer.start_line = False
        return t

    def t_pre_LINK(self, t):
        r'\[\s*\S+\s*]'
        self.lexer.start_line = False
        return t

    def t_pre_ORPHAN(self, t):
        r'((http://|ftp://)\S+)'
        self.lexer.line_start = False
        if not self.lexer.code_block:
            url = t.value.replace(",","")
            t.value = "<%s>" % url
        self.lexer.start_line = False
        return t
    
    def t_ANY_CODE_OPEN(self, t):
        r'<code>'
        self.lexer.code_block = True
        return t
    
    def t_ANY_CODE_CLOSE(self, t):
        r'</code>'
        self.lexer.code_block = False
        return t
    
    def t_ANY_TAG_OPEN(self, t):
        r'<\w+>'
        self.lexer.start_line = False
        return t

    def t_ANY_TAG_CLOSE(self, t):
        r'</\w+>'
        self.lexer.start_line = False
        return t

    def t_post_COMMAND(self, t):
        r'\\\w+\s+[a-zA-Z0-9_,\-]*'
        if not self.lexer.code_block:
            t.value = parse_command(t.value)
        self.lexer.start_line = False
        return t

    def t_ANY_TEXT(self, t):
        r'(\S+)'
        self.lexer.line_start = False
        return t

    # A string containing ignored characters (spaces and tabs)
    t_ANY_ignore  = '\t'

    # Error handling rule
    def t_ANY_error(self, t):
        print "Illegal character '%s'" % t.value[0]
        t.lexer.skip(1)
         
    # Build the lexer
    def build(self,**kwargs):
        self.lexer = lex.lex(module=self, **kwargs)
        self.lexer.code_block = False
        self.lexer.start_line = True
    
    # Test it output
    def test(self,data):
        self.lexer.input(data)
        while True:
             tok = self.lexer.token()
             if not tok: 
                 break
             print tok,  self.lexer.code_block, tok.value

def user_html(vals):
    from main.server import models
    posts = models.User.objects.filter(id__in=vals).select_related('user__profile')
    patt  = '<a href="%s">%s</a>'
    coll  = [patt % (u.profile.get_absolute_url(), u.profile.display_name) for u in posts ]
    return ", " .join(coll)
    
def tag_html(vals):
    tags = "+".join(vals)
    patt  = '<a href="/show/tag/%s/">%s</a>' %(tags, tags)
    return patt
    
def post_html(vals):
    from main.server import models
    posts = models.Post.objects.filter(id__in=vals)
    patt  = '<a href="%s">%s</a>'
    coll  = [patt % (p.get_absolute_url(), p.title) for p in posts ]
    return ", " .join(coll)

def gist_html(vals):
    patt  = '<div><script src="http://gist.github.com/%s.js"></script></div>'
    coll  = [patt % v for v in vals ]
    return ", " .join(coll)

def search_html(vals):
    html = """
<div id="cse" style="width: 100%;">Loading</div>
<script src="http://www.google.com/jsapi" type="text/javascript"></script>
<script type="text/javascript"> 
  google.load('search', '1', {language : 'en', style : google.loader.themes.GREENSKY});
  google.setOnLoadCallback(function() {
    var customSearchOptions = {};  var customSearchControl = new google.search.CustomSearchControl(
      '003596843386727440968:raditaczxza', customSearchOptions);
    customSearchControl.setResultSetSize(google.search.Search.FILTERED_CSE_RESULTSET);
    customSearchControl.draw('cse');
  }, true);
</script>
</div>
"""
    return html
    
def youtube_html(vals):
    patt = r'''
    <div>
        <iframe width="560" height="315" src="http://www.youtube.com/embed/%s" frameborder="0" allowfullscreen></iframe>
    </div>
    '''
    coll  = [patt % v for v in vals ]
    return ", " .join(coll)

CMD_MAP = {
    '\\post':post_html,
    '\\user':user_html,
    '\\gist':gist_html,
    '\\tag':tag_html,
    '\\tags':tag_html,
    '\\search': search_html,
    '\\youtube':youtube_html,
}

def parse_command(text):
    "parses a shortcut"
    cmd, vals = text.split()
    vals = vals.split(",")[:10] # sanity check
    out = CMD_MAP.get(cmd, lambda x:text)(vals)
    return out

def process(text, state='post'):
    "Process text"
    try:
        p = MainLexer()
        p.build()  
        p.lexer.input(text)
        p.lexer.begin(state)
        c = []
        while 1:
            tok = lex.token()
            if not tok: break
            c.append(tok.value)
            #print tok,  p.lexer.code_block, tok.value
        text = "".join(c)
    except Exception, exc:
        print "*** textparser error: %s" % exc
    return text

if __name__ == '__main__':
    # Test it out
    text = '''Dear all,

http://www.biostars1.org

I am using MACS and to analyze chip-seq data. http://www.biostars2.org \gist 2059 \post 10 The PeakSplitter
here is something. http://www.biostars3.org ftp://something?a=b&b=c:5
    http://www.biostars5.org
    ok this is still inside a codeblock http://www.biostars6.org
[  http://www.biostars7.org      ](ok)
  http://www.biostars8.org not indented sufficiently
[something ]:   http://www.biostars.org
ok this is outside a codeblock http://www.nd.edu
\\gist 123,456
\\youtube 1abcde
<p>\\post 3,4,5</p>
'''

    m = MainLexer()
    m.build()  
    #m.test(text)
   
    out = process(text, state='pre')
    print out

    print '-' * 30

    out = process(text, state='post')

    #print out

   