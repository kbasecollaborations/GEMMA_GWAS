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

var bonferroni = -1*Math.log10((0.05/ind));


d3.tsv(inputs[0], function(error, data) {
  data.forEach(function(d) {
    d.bp = +d.BP;
    d["chr"] = d["CHR"];
    d.id = +d.POS;
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

  var color = d3.scale.category10();

  var tip = d3.tip()
      .attr("class", "d3-tip")
      .offset([-10, 0])
      .html(function(d) {
        tooltiphtml = "Chromosome : "+ d["chr"] + "<br/>";
        tooltiphtml += "Chr Pos:" + d.bp + "<br/>";
        tooltiphtml += yCat + ": " + d.logp;

        if(d.gene != NaN && d.gene != "NA") {
          tooltiphtml += '<br/><span style="color: #00ff2a;">GENE: '+d.gene+'</span>';
        } else if(d.neighborgene != NaN && d.neighborgene != "NA") {
          tooltiphtml += '<br/><span style="color: #00ff2a;">NEAREST GENE: '+d.neighborgene+'</span>';
        }

        if (d.genefunction != NaN && d.genefunction != "NA") {
          tooltiphtml += '<p style="width: 300px;margin:0;padding:0;color: #ff0061;word-wrap: break-word;">GENE FUNCTION:<br/>'+d.genefunction+'</p>';
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
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


  svg.append('line')
    .attr('x1', x(0))
    .attr('x2', x(xMax))
    .attr('y1', y(bonferroni))
    .attr('y2', y(bonferroni))
    .attr('stroke', '#00b500');

  svg.append('text')
    .attr("x", x(0))
    .attr("y", y(bonferroni))
    .text("Bonferroni")
    .attr("font-family", "sans-serif")
    .attr("font-size", "24px")
    .attr("fill", "#00b500");

  svg.call(tip)
     .call(zoomBeh);

  svg.append("rect")
      .attr("width", width)
      .attr("height", height)
      .attr("stroke", "#00b500");

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

  svg.append("text")
        .attr("x", (width / 2))
        .attr("y", 0 - (margin.top / 2))
        .attr("text-anchor", "middle")
        .style("font-size", "16px")
        .style("text-decoration", "underline")
        .text("Manhattan Plot");

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

