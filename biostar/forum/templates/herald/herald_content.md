{% load forum_tags %}
{% load humanize %}
{% load static %}

The **Biostar Herald** publishes user submitted links of bioinformatics relevance. It aims to provide a summary of interesting and relevant information you may have missed. You too can submit [links here](/herald/).

This edition of the Herald was brought to you by contribution from {% for author in authors %} [{{author.profile.name}}]({{base_url}}{{ author.profile.get_absolute_url }}),
{% endfor %} and was edited by {% for editor in editors %} [{{editor.profile.name}}]({{base_url}}{{ editor.profile.get_absolute_url }}),
{% endfor %}

{% for herald in heralds %}

---

{% if herald.title %}
### [{{herald.title}}]({{herald.url}}) ({{herald.domain}})
{% else %}
### {{herald.url}}
{% endif %}

{% if herald.text %}
{{herald.text|safe}}
{% endif %}

submitted by: [{{herald.author.profile.name}}]({{base_url}}{{ herald.author.profile.get_absolute_url }})

{% endfor %}

---

Want to get the **Biostar Herald** in your email? Who wouldn't? Sign up righ'ere: <a class="herald-sub">toggle subscription</a>






