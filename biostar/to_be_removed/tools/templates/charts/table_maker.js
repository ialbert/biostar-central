
    google.charts.load("current", {packages:["table"]});
    google.charts.setOnLoadCallback({{params.table_id}});

    function {{params.table_id}}() {

        var data = new google.visualization.DataTable();
        {% for  name, value in params.columns %}
        data.addColumn('{{name}}', '{{value}}');
        {% endfor %}
        data.addRows([
        {%  for row in params.rows %}
        [ {{row |safeseq|join:", " }} ],
        {%  endfor %}
        ]);

        var options = {
        {{params.options|safe}}
        };

    var table = new google.visualization.{{params.type}}(document.getElementById('{{params.table_id}}'));
    table.draw(data, options);

    }
