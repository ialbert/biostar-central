from django.core.management.base import BaseCommand
from whoosh.fields import Schema, TEXT, KEYWORD, ID, STORED
from whoosh.analysis import StemmingAnalyzer




class Command(BaseCommand):


    help = 'Batch update the search index with current database'

    def add_arguments(self, parser):
        parser.add_argument('--root', required=True, help="root")

    def handle(self, *args, **options):

        schema = Schema(from_addr=ID(stored=True),
                        to_addr=ID(stored=True),
                        subject=TEXT(stored=True),
                        body=TEXT(analyzer=StemmingAnalyzer()),
                        tags=KEYWORD)
