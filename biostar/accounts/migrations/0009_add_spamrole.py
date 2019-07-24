
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('accounts', '0008_message_body'),
    ]

    operations = [
        migrations.AlterField(
            model_name='profile',
            name='role',
            field=models.IntegerField(choices=[(0, 'Reader'), (1, 'Moderator'), (2, 'Admin'), (3, 'Blog User'), (4, 'Spammer')], default=0),
        ),
    ]
