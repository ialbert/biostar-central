from biostar.tools import render
import sys, csv

# TODO :  Move template to a different file.

template = '''
<html>
  <head>
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">

      // Load Charts and the corechart, barchart and map packages.
      google.charts.load('current', {'packages':['corechart','map']});

      // Draw the pie chart and bar chart with the same data.
      google.charts.setOnLoadCallback(drawChart);
      
      // Draw map.
      google.charts.setOnLoadCallback(drawMap);

      function drawChart() {

        var data = new google.visualization.DataTable();
        data.addColumn('string', 'Species');
        data.addColumn('number', 'Counts');
        data.addRows([
         {% for name, value in items %}
            ['{{name}}' ,{{value}}],
          {% endfor %}
        ]);

        var piechart_options = {title:'Unique reads in species',
                       chartArea:{left:20,top:0,width:'80%',height:'75%'},
                       legend:{position:'left'},
                       };
        var piechart = new google.visualization.PieChart(document.getElementById('piechart_div'));
        piechart.draw(data, piechart_options);

        var barchart_options = {title:'Unique reads in species',
                       width :800,
                       height:500,
                       legend: 'none'};
        var barchart = new google.visualization.BarChart(document.getElementById('barchart_div'));
        barchart.draw(data, barchart_options);
      }
      
      function drawMap() {
      var data = google.visualization.arrayToDataTable([
        ['Lat', 'Long', 'Name'],
        [41.015908, -77.53246, 'Lamar PA']
      ]);

    var options = {
      zoomLevel: 6,
      showTooltip: true,
      showInfoWindow: true
    };

    var map = new google.visualization.Map(document.getElementById('mapchart_div'));

    map.draw(data, options);
  };
                
      
</script>
<body>
    <h2> Classification Results<h2>
    
    <div id="barchart_div" style="width: 900px; height: 500px;"></div> 
    <div id="piechart_div" style="width: 900px; height: 500px;"></div>    
     <div id="mapchart_div" style="width: 900px; height: 500px;"></div>
    
  </body>
</html>
'''


#fname = "report.txt"
# template_file = "classify_results.html"
# template = open(template_file).read()


fname = sys.argv[1]
items = []

# Read centrifuge report file.
with open(fname) as csvfile:
    reader = csv.DictReader(csvfile, delimiter="\t")
    for row in reader:
        species = row['name']
        reads = row['numUniqueReads']
        items.append((species, reads))

# Sort by counts.
# items.sort(key=lambda x: int(x[1]))

# Sort by species name.
items.sort(key=lambda x: x[0])

# render html
context = dict(items=items)
html = render.render_data(template_txt=template, context=context)
print(html)
