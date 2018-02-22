
    google.charts.load("current", {packages:["corechart"]});
    google.charts.setOnLoadCallback({{params.chart_id}});

    function {{params.chart_id}}() {
        var data = google.visualization.arrayToDataTable([

            ['{{ params.xlabel }}', '{{ params.ylabel }}' ],
        {%  for name, value in params.data %}
            [ '{{name|safe}}', {{ value|safe }} ],
        {%  endfor %}

        ]);

    var options = {
        {{params.options|safe}}
        };

    var chart = new google.visualization.{{params.type}}(document.getElementById('{{params.chart_id}}'));
    chart.draw(data, options);

}
