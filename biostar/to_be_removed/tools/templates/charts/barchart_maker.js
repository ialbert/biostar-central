
    google.charts.load("current", {packages:["corechart"]});
    google.charts.setOnLoadCallback({{params.chart_id}});

    function {{params.chart_id}}() {
        var data = google.visualization.arrayToDataTable([

            [ {{params.header|safeseq|join:", " }} ],
        {%  for row in params.data %}
            [ {{row |safeseq|join:", " }} ],
        {%  endfor %}

        ]);

    var options = {
        {{params.options|safe}}
        };

    var chart = new google.visualization.{{params.type}}(document.getElementById('{{params.chart_id}}'));
    chart.draw(data, options);

}
