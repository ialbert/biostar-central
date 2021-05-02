Here are the latest things shared amongst our user base.

{% for herald in heralds %}
  
   - {{herald.user.profile.name}} shared {{herald.url}}
   
{% endfor %}