{% load forum_tags %}
{% load humanize %}
{% load static %}


{% for herald in heralds %}

 {{herald.url}} shared by {{base_url}}{{ herald.author.profile.get_absolute_url }}
 
{% endfor %}
