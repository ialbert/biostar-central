import logging

from django.core.management.base import BaseCommand
from biostar.forum import markdown

logger = logging.getLogger('engine')

test = """

<b>foo</b>

1 > 0

    1 > 0 
 
"""



test2 = """

<b> foo2 </b>

<script> bar() </script>

"""

gist_test = "https://gist.github.com/afrendeiro/6732a46b949e864d6803"

class Command(BaseCommand):
    help = 'Used to test markdown rendering'

    def handle(self, *args, **options):
        #import markdown2

        #html1 = markdown2.markdown(test)
        #print(html1)

        #print('MISTUNE', '-'*50)
        html = markdown.parse(test, escape=False, clean=False, allow_rewrite=False)
        print(html)
