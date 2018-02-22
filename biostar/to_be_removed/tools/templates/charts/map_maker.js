
    google.charts.load("current", {packages:["map"]});
    google.charts.setOnLoadCallback({{params.chart_id}});

    function {{params.chart_id}}() {
        var data = google.visualization.arrayToDataTable([

            ['Lat', 'Long', 'Name' ],
        {%  for latitude, longitude, name in params.data %}
            [ {{latitude|safe}}, {{ longitude|safe }}, '{{name|safe}}' ],
        {%  endfor %}

        ]);

    var options = {
        {{params.options|safe}}
        };

    var chart = new google.visualization.{{params.type}}(document.getElementById('{{params.chart_id}}'));
    chart.draw(data, options);

}
