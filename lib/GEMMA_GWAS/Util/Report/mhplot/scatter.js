var margin = { top: 50, right: 300, bottom: 50, left: 50 },
    outerWidth = 1050,
    outerHeight = 500,
    width = outerWidth - margin.left - margin.right,
    height = outerHeight - margin.top - margin.bottom;
    width = 960 - margin.left - margin.right,
    height = 400 - margin.top - margin.bottom;

var x = d3.scale.linear()
    .range([0, width]).nice();

var y = d3.scale.linear()
    .range([height, 0]).nice();

var xCat = "POS",
    yCat = "P",
    rCat = "CHR",
    colorCat = "CHR";

/*
d3.csv("cereal.csv", function(data) {
  data.forEach(function(d) {
    d.Calories = +d.Calories;
    d.Carbs = +d.Carbs;
    d["Cups per Serving"] = +d["Cups per Serving"];
    d["Dietary Fiber"] = +d["Dietary Fiber"];
    d["Display Shelf"] = +d["Display Shelf"];
    d.Fat = +d.Fat;
    d.Potassium = +d.Potassium;
    d["Protein (g)"] = +d["Protein (g)"];
    d["Serving Size Weight"] = +d["Serving Size Weight"];
    d.Sodium = +d.Sodium;
    d.Sugars = +d.Sugars;
    d["Vitamins and Minerals"] = +d["Vitamins and Minerals"];
  });
*/


d3.tsv(inputs[0], function(error, data) {
  data.forEach(function(d) {
    d.bp = +d.BP;
    d["chr"] = d["CHR"];
    // alert(d.CHR);
    // alert(d.chr);
    d.id = +d.POS;
    //alert(d.id);
    d.p = +d.P;
    d.logp = -Math.log10(d.p);
    d["gene"] = d["GENEID"];
    d["neighborgene"] = d["NEIGHBORGENE"];
    d["genefunction"] = d["FUNCTION"];
  });


  var xMax = d3.max(data, function(d) { return d.id; }) * 1.05,
      xMin = d3.min(data, function(d) { return d.id; }),
      xMin = xMin > 0 ? 0 : xMin,
      yMax = d3.max(data, function(d) { return d.logp; }) * 1.05,
      yMin = d3.min(data, function(d) { return d.logp; }),
      yMin = yMin > 0 ? 0 : yMin;


  x.domain([xMin, xMax]);
  y.domain([yMin, yMax]);

  var xAxis = d3.svg.axis()
      .scale(x)
      .orient("bottom")
      .tickSize(-height);

  var yAxis = d3.svg.axis()
      .scale(y)
      .orient("left")
      .tickSize(-width);
  //xAxis.tickValues([10000000, 20000000, 30000000]);
  //.tickValues(d3.range(1000000, 2000000, 4));
  var color = d3.scale.category10();

  var tip = d3.tip()
      .attr("class", "d3-tip")
      .offset([-10, 0])
      .html(function(d) {
        tooltiphtml = "Chromosome : "+ d["chr"] + "<br/>";
        tooltiphtml += "Chr Pos:" + d.bp + "<br/>";
        tooltiphtml += yCat + ": " + d.logp;

        if(d.gene != NaN && d.gene != "NA") {
          tooltiphtml += '<br/><span style="color: #00ff2a;">GENE: '+d.gene+'</
        } else if(d.neighborgene != NaN && d.neighborgene != "NA") {
          tooltiphtml += '<br/><span style="color: #00ff2a;">NEAREST GENE: '+d.
        }

        if (d.genefunction != NaN && d.genefunction != "NA") {
          tooltiphtml += '<p style="width: 300px;margin:0;padding:0;color: #ff0
        }

        return tooltiphtml
      });

  var zoomBeh = d3.behavior.zoom()
      .x(x)
      .y(y)
      .scaleExtent([0, 500])
      .on("zoom", zoom);

  var svg = d3.select("#scatter")
    .append("svg")
      .attr("width", outerWidth)
      .attr("height", outerHeight)
    .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
      .call(zoomBeh);

  svg.call(tip);

  svg.append("rect")
      .attr("width", width)
      .attr("height", height);

/*
  svg.append("g")
      .classed("x axis", true)
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
    .append("text")
      .classed("label", true)
      .attr("x", width)
      .attr("y", margin.bottom - 10)
      .style("text-anchor", "end")
      .text("Chromosomes");


  svg.append("g")
	        .attr("class", "x axis")
	        .attr("transform", "translate(0," + height + ")")
	        .call(xAxis)
	        .selectAll("text")
	            .style("text-anchor", "end")
	            .attr("dx", "-.8em")
	            .attr("dy", ".15em")
		    .attr("transform", "rotate(-90)" );


svg.append("g")
	                .attr("class", "x axis")
	                .attr("transform", "translate(0," + height + ")")
	                .call(xAxis)
	                .selectAll("text")
	                 .style("text-anchor", "end")
	                    .attr("dx", "-.8em")
	                    .attr("dy", ".15em")
	                    .attr("transform", "rotate(-90)" )
*/
  svg.append("g")
      .classed("y axis", true)
      .call(yAxis)
    .append("text")
      .classed("label", true)
      .attr("transform", "rotate(-90)")
      .attr("y", -margin.left)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text("-log10(p-value)");

  var objects = svg.append("svg")
      .classed("objects", true)
      .attr("width", width)
      .attr("height", height);

  objects.append("svg:line")
      .classed("axisLine hAxisLine", true)
      .attr("x1", 0)
      .attr("y1", 0)
      .attr("x2", width)
      .attr("y2", 0)
      .attr("transform", "translate(0," + height + ")");

  objects.append("svg:line")
      .classed("axisLine vAxisLine", true)
      .attr("x1", 0)
      .attr("y1", 0)
      .attr("x2", 0)
      .attr("y2", height);

  objects.selectAll(".dot")
       .data(data)
       .enter().append("circle")
       .classed("dot", true)
       .attr("r", 3.5)
       .attr("transform", transform)
       .style("fill", function(d) { return color(d[colorCat]); })
       .on("mouseover", tip.show)
       .on("mouseout", tip.hide);

/*
  objects.selectAll(".dot")
      .data(data)
    .enter().append("circle")
      .classed("dot", true)
      .attr("r", function (d) { return 6 * Math.sqrt(d[rCat] / Math.PI); })
      .attr("transform", transform)
      .style("fill", function(d) { return color(d[colorCat]); })
      .on("mouseover", tip.show)
      .on("mouseout", tip.hide);
*/
  svg.append("text")
        .attr("x", (width / 2))
        .attr("y", 0 - (margin.top / 2))
        .attr("text-anchor", "middle")
        .style("font-size", "16px")
        .style("text-decoration", "underline")
        .text("Manhattan Plot");
/*
  var legend = svg.selectAll(".legend")
      .data(color.domain())
    .enter().append("g")
      .classed("legend", true)
      .attr("transform", function(d, i) { return "translate(0," + i * 20 + ")"; });

  legend.append("circle")
      .attr("r", 3.5)
      .attr("cx", width + 20)
      .attr("fill", color);

  legend.append("text")
      .attr("x", width + 26)
      .attr("dy", ".35em")
      .text(function(d) { return  d; });
*/

  function zoom() {
    svg.select(".x.axis").call(xAxis);
    svg.select(".y.axis").call(yAxis);

    svg.selectAll(".dot")
        .attr("transform", transform);
  }

  function transform(d) {
    return "translate(" + x(d[xCat]) + "," + y(d.logp) + ")";
  }
});

