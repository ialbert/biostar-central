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
       'TAG_OPEN',
       'TAG_CLOSE',
       'COMMAND',
       'TEXT',
    )

    def t_ANY_NEWLINE(self, t):
        r'\n|\r\n'
        self.lexer.code_block = False
        self.lexer.start_line = True
        return t

    def t_ANY_SPACE(self, t):
        r'(\s+)'
        if len(t.value) == 4 and self.lexer.start_line:
            self.lexer.code_block = True
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
            t.value = "<%s>" % t.value
        self.lexer.start_line = False
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
    '\\gist':gist_html,
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

http://www.nd.edu

I am using MACS and to analyze chip-seq data. http://www.nd.edu \gist 2059 \post 10 The PeakSplitter
here is something. http://www.biostars.org ftp://something?a=b&b=c:5
    http://www.nd.edu
    ok this is still inside a codeblock http://www.nd.edu
[  http://www.biostars.org      ](ok)
  http://www.x.com
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

    print out

   