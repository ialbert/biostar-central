from rest_framework import serializers

from .models import Vote, Post


class VoteSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = Vote
        fields = ('id', 'url', 'date', 'author', 'post')
                #TODO type


class PostSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = Post
        fields = ('id', 'url', 'title', 'author', 'votes')
                #TODO add other fields