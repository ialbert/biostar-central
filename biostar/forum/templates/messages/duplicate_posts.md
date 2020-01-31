
Hello {{ user.profile.name }}!

Posts similar to yours can already be found at:

{% for dupe_url in dupes %}
- [{{ dupe_url }}]({{dupe_url}})
{% endfor %}

We encourage you to look at the sources listed above for more answers. 

If you disagree with this please tell us why in a reply below. We'll be happy to talk about it.

Cheers!

{%  if comment %}
 PS: {{ comment }}
{%  endif %}