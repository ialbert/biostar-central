
**Congratulations!** You have won the **[{{ award.badge.name }}]({% url 'badge_view' award.badge.uid %})** <i class="{{ award.badge.icon }}"></i>
badge for the notable accomplishment of: *{{ award.badge.desc }}*.{% if award.post and award.post.root %}
You have earned this honor for: [{{ award.post.title}}]({{ award.post.get_absolute_url }})
{% endif %}

