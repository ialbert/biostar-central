{% load forum_tags %}
{% load humanize %}
{% load static %}

{% include 'herald/herald_top.md' %}

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

{% include 'herald/herald_bottom.md' %}