{% load engine_tags %}

{% spaceless %}
{% if data.items %}
Parameters used during the run:
{% for key, obj in data.items %}{% if obj.label %}- {{ obj.label }}:  {% if obj.name %}__{{ obj.name }}__{% else %}__{{ obj.value }}__{% endif %}{% endif %}
{% endfor %}
{% endif %}
{% endspaceless %}
