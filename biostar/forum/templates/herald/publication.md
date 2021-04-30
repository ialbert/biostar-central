
{% for herald in heralds %}
   
   - {{herald.user.profile.name}} suggest {{herald.url}}
   
{% endfor %}