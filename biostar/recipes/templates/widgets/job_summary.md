{% load engine_tags %}
#### Parameters used during the run:
{% if data.items %}{% for key, obj in data.items %}{% if obj.label %} - {{ obj.label }}:{% if obj.name %}__{{ obj.name }}__{% else %}__{{ obj.value }}__{% endif %}
{% endif %}{% endfor %}
{% else %}No user parameters were set.{% endif %}
