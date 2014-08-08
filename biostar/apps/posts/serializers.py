from rest_framework import serializers

from .models import Vote


class VoteSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = Vote
        fields = ('id', 'url', 'date')
                #author
                #post
                #type