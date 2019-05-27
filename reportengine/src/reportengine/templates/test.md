This is a test with title {{title}}
===================================

{{ resource }}

{% for i in range(5) %}

{{ otherresource(param=i) }}

{% endfor %}

{{ otherresource }}
{{ otherresource (param=3) }}

