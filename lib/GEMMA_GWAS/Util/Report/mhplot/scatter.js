// Create data
/*
function randomData(samples) {
    var data = [],
        random = d3.randomNormal();
        
    for (i = 0; i < samples; i++) {
        data.push({
            x: random(),
            y: random()
        });
    }
    return data;
}

var data = randomData(300);
*/
// tooltip mouseover event handler

d3.tsv(inputs[0], function(error, data) {
    data.forEach(function(d) {
        d.bp = d["POS"];
        d["chr"] = d["CHR"];    
        d.p = d["P"];
        d["gene"] = d["GENEID"];
        d["neighborgene"] = d["NEIGHBORGENE"];
        d["genefunction"] = d["FUNCTION"];
        d.x = Number(+d["BP"]);
        d.y = (-1*Math.log10(d.p));
    });

    var colorCat = "CHR";
    var color = d3.scaleOrdinal(d3.schemeCategory20c);

    var margin = { top: 20, right: 20, bottom: 30, left: 30 };
    width = 900 - margin.left - margin.right,
    height = 480 - margin.top - margin.bottom;

    var tooltip = d3.select("body").append("div")
      .attr("id", "tooltip")

    var x = d3.scaleLinear()          
          .range([0, width])
          .nice();

    var y = d3.scaleLinear()
        .range([height, 0]);

    var xAxis = d3.axisBottom(x).ticks(12),
        yAxis = d3.axisLeft(y).ticks(12 * height / width);

    var brush = d3.brush().extent([[0, 0], [width, height]]).on("end", brushended),
        idleTimeout,
        idleDelay = 350;

    var svg = d3.select("body").append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("g")
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var clip = svg.append("defs").append("svg:clipPath")
        .attr("id", "clip")
        .append("svg:rect")
        .attr("width", width )
        .attr("height", height )
        .attr("x", 0) 
        .attr("y", 0); 

    var xExtent = d3.extent(data, function (d) { return d.x; });
    var yExtent = d3.extent(data, function (d) { return d.y; });
    x.domain(d3.extent(data, function (d) { return d.x; })).nice();
    y.domain(d3.extent(data, function (d) { return d.y; })).nice();

    var scatter = svg.append("g")
         .attr("id", "scatterplot")
         .attr("clip-path", "url(#clip)");

    scatter.append("g")
        .attr("class", "brush")
        .call(brush);
        
    scatter.selectAll(".dot")
        .data(data)
        .enter()
        .append("circle")
        .attr("class", "dot")
        .attr("r", 4)
        .attr("cx", function (d) { return x(d.x); })
        .attr("cy", function (d) { return y(d.y); })
        .attr("opacity", 0.5)
        .style("fill", function (d) { return color(d[colorCat]); })
        .style("pointer-events","visible")
        .on('mouseover', d => {
           tooltip.transition()
             .duration(100)
             .style('opacity', .9);
           tooltip.html(function() { 
                tooltiphtml = "<p>Chromosome: "+ d["chr"] + "</p>";
                tooltiphtml += "<p>Position: " + d.bp + "</p>";
                tooltiphtml += "<p>Pvalue: " + d.y.toFixed(2) + "</p>";

                if(d.gene != NaN && d.gene != "NA" && d.gene != undefined) {
                  tooltiphtml += '<p>GENE: '+d.gene+'</p>';
                } else if(d.neighborgene != NaN && d.neighborgene != "NA" && d.neighborgene != undefined) {
                  tooltiphtml += '<p>NEAREST GENE: '+d.neighborgene+'</p>';
                }

                if (d.genefunction != NaN && d.genefunction != "NA" && d.genefunction != undefined) {
                  tooltiphtml += '<p>GENE FUNCTION: '+d.genefunction+'</p>';
                }

                return tooltiphtml
           })
             .style('left', `${d3.event.pageX + 2}px`)
             .style('top', `${d3.event.pageY - 18}px`);
        })
        .on('mouseout', () => {
           tooltip.transition()
             .duration(400)
             .style('opacity', 0);
        });

    // x axis
    svg.append("g")
       .attr("class", "x axis")
       .attr('id', "axis--x")
       .attr("transform", "translate(0," + height + ")")
       .call(xAxis);

    svg.append("text")
     .style("text-anchor", "end")
        .attr("x", width)
        .attr("y", height - 8)
     .text("Chromosome")
        .attr('class', "stroke-axis-text");

    // y axis
    svg.append("g")
        .attr("class", "y axis")
        .attr('id', "axis--y")
        .call(yAxis);

    svg.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", "1em")
        .style("text-anchor", "end")
        .text("-log10(p)")
            .attr('class', "stroke-axis-text");

    function brushended() {
        var s = d3.event.selection;
        if (!s) {
            if (!idleTimeout) return idleTimeout = setTimeout(idled, idleDelay);
            x.domain(d3.extent(data, function (d) { return d.x; })).nice();
            y.domain(d3.extent(data, function (d) { return d.y; })).nice();
        } else {
            
            x.domain([s[0][0], s[1][0]].map(x.invert, x));
            y.domain([s[1][1], s[0][1]].map(y.invert, y));
            scatter.select(".brush").call(brush.move, null);
        }
        zoom();
    }

    function idled() {
        idleTimeout = null;
    }

    function zoom() {
        var t = scatter.transition().duration(750);
        svg.select("#axis--x").transition(t).call(xAxis);
        svg.select("#axis--y").transition(t).call(yAxis);
        scatter.selectAll("circle").transition(t)
        .attr("cx", function (d) { return x(d.x); })
        .attr("cy", function (d) { return y(d.y); });
    }
});