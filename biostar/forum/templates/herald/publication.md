{% load forum_tags %}

Here are the latest things shared by our user base.

{% for herald in heralds %}
  
   - <img class="ui avatar image" src="{% gravatar user=herald.user size=50 %}"/> {{herald.user.profile.name}} shared {{herald.url}}
   
{% endfor %}