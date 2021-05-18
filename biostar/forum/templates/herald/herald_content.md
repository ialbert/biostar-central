{% load forum_tags %}
{% load humanize %}
{% load static %}

The **Biostar Herald** publishes user submitted links with bioinformatics relevance. It attempts to be as a snapshot of interesting and relevant information you may have missed, served in a "digest" format. You too can submit [interesting links here](/herald/) and with that you can help publish the next issue.

This edition of the Herald was brought to you by contribution from *TODO (authors)* and was edited by *TODO (mods)*

{% for herald in heralds %}

<div class="ui divider"></div>

[{{herald.url|truncatechars:100}}]({{herald.url}})

{% if herald.text %}
{{herald.text}}
{% endif %}

submitted by: [{{herald.author.profile.name}}]({{base_url}}{{ herald.author.profile.get_absolute_url }})

{% endfor %}

<div class="ui divider"></div>

Want to get the **Biostar Herald** in your email? Who wouldn't? Sign up righ'ere: *link goes here*






