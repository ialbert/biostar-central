template = '''<html>
  <head>\
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load("current", {packages:["corechart"]});
      google.charts.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = google.visualization.arrayToDataTable([
          {% for name, value in items %}
          {% endfor %}

        var options = {
          title: 'Lengths of dinosaurs, in meters',
          legend: { position: 'none' },
        };

        var chart = new google.visualization.Histogram(document.getElementById('chart_div'));
        chart.draw(data, options);
      }
    </script>
  </head>
  <body>
    <h2> Classification Results<h2>
    
    <div id="chart_div" style="width: 900px; height: 500px;"></div>
  </body>
</html>
'''

items = get_template()

contxext = dict(items=items)

html = render(tempalte, context)

print (html)
