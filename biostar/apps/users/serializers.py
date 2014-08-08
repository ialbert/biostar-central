from datetime import datetime
from rest_framework import serializers

from biostar.apps.posts.models import Vote
from .models import User


class UserSerializer(serializers.HyperlinkedModelSerializer):
    date_joined = serializers.DateTimeField(source='profile.date_joined')
    joined_days_ago = serializers.SerializerMethodField('get_joined_days_ago')
    vote_count = serializers.SerializerMethodField('get_vote_count')
    #questions = ...
    #answers = ...

    class Meta:
        model = User
        fields = ('id', 'url', 'email', 'name', 'last_login', 'date_joined', 'joined_days_ago',
                  'vote_count')

    def get_joined_days_ago(self, obj):
        return (datetime.now().date() - obj.profile.date_joined.date()).days

    def get_vote_count(selfself, obj):
        return Vote.objects.filter(author=obj).count()