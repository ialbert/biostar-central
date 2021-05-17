{% load forum_tags %}
{% load humanize %}
{% load static %}

The **Biostar Herald** publishes user submitted links with bioinformatics relevance. [Submit your links here.](/herald/)

This edition of the Herald was brought to you by: *name of the peeps sharing*

Edited by:  *name of the mods accepting links*

{% for herald in heralds %}

<div class="ui divider"></div>

Submitted by: [{{herald.author.profile.name}}]({{base_url}}{{ herald.author.profile.get_absolute_url }})

{% if herald.text %}
{{herald.text}}
{% endif %}

* [{{herald.url|truncatechars:50}}]({{herald.url}})

{% endfor %}

<div class="ui divider"></div>

Want to get the **Biostar Herald** in your email? Who wouldn't? Sign up righ'ere: *link goes here*






