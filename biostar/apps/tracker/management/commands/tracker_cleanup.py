# -*- encoding: utf-8 -*-
from __future__ import (unicode_literals, division, absolute_import)

from django.core.management.base import BaseCommand, CommandError

class Command(BaseCommand):
    args = ''
    help = 'Cleanup peer tables'

    def handle(self, *args, **options):
        self.stdout.write('Cleaning up dead peers\n')
        cleanup_peers()
        self.stdout.write('Cleaning up orphaned peerinfo records\n')
        cleanup_peerinfo()
        self.stdout.write('Done\n')
