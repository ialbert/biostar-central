import os
from django.db.models.signals import post_migrate
from django.core.files import File
from django.apps import AppConfig
from django.conf import settings


class EngineConfig(AppConfig):
    name = 'biostar.recipes'

    def ready(self):
        # Triggered upon app initialization.
        post_migrate.connect(init_commands, sender=self)
        pass


def init_commands(sender, **kwargs):
    from biostar.recipes.models import CommandType, Command
    from biostar.accounts.models import User

    image_root = os.path.join(settings.BASE_DIR, 'initial')

    inital_cmds = (('bash', 'bashicon.png', dict(label='Print all executed commands.', command='set -uexo pipefail')),
                   ('bash', '', dict(label='Print hello.', command='echo HELLO', image='bashicon.png')),
                   ('r', 'r_icon.png', dict(label='Generate a sample.', command='sample(1:3, size=1000, prob=c(.30,.60,.10))')),
                   ('matlab', 'matlab.png', dict(label=' Open a file in read mode.', command='fopen(filename,'r')')),
                   ('kraken', 'kraken.png', dict(label=' Add file to the database.', command='kraken2-build --add-to-library fish-accession.fa -db db 2')),
                   ('centrifuge', 'centrifuge.png', dict(label='Build the index.', command='centrifuge-build -p $N --conversion-table $TABLE --taxonomy-tree $NODES  --name-table $NAMES  $REFERENCE $INDEX')),
                   ('qiime', 'qiime2.png', dict(label='Convert taxonomy file to qiime 2 artifact.', command=" qiime tools import --input-path $REF_FASTA --output-path $REFERENCE --type 'FeatureData[Sequence]'"))
                   )
    owner = User.objects.filter(is_superuser=True).first()
    for name, image, data in inital_cmds:

        cmd_type = CommandType.objects.filter(name=name).first()

        if cmd_type is None:
            cmd_type = CommandType.objects.create(name=name)
            if image:
                image = File(open(os.path.join(image_root, image), 'rb'))
                cmd_type.image = image
                cmd_type.save()

        help_text, command = data['label'], data['command']

        Command.objects.create(help_text=help_text, command=command,
                               type=cmd_type, owner=owner)

