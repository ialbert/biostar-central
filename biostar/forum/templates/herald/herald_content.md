{% load forum_tags %}
{% load humanize %}
{% load static %}

<div class="ui basic segment">

<h2 class="ui header">
<img src="{% static 'images/news-herald.png' %}">
<div class="content">
Biostar Herald
<div class="sub header">Share bioinformatics resources from across the web.</div>
</div>
</h2>

</div>

{% for herald in heralds %}
 
<div class="ui container">

<img class="ui avatar image" src="{% gravatar user=herald.author size=50 %}"/>
<a href="{{ herald.author.profile.get_absolute_url }}">
{{herald.author.profile.name}}
</a>&bull; 
<em>{{ herald.date|naturaltime }}</em> &bull;
accepted by <a href="{{ herald.editor.profile.get_absolute_url }}">
 {{herald.editor.profile.name}}
</a>

<div class='ui basic segment'>{{herald.text}}</div>
</div>

{% endfor %}