var StreamHist = require('./index.js').StreamHist,
    rand = require('randgen'),  // Just for some random data
    d3 = require("d3"),
    doc = require("jsdom").jsdom();

// Create random Gaussian data and push it into a StreamHist object
var data = rand.rvnorm(10000)
var hist = StreamHist(20, false, 0, 5000).push(data);

// Employ conventional margins (http://bl.ocks.org/3019563)
var margin = {top: 20, right: 20, bottom: 50, left: 40},
    width = 300 - margin.left - margin.right,
    height = 300 - margin.top - margin.bottom;

var x = d3.scaleLinear()
    .range([0, width])
    .domain(hist.limits());

var y = d3.scaleLinear()
    .range([height, 0])
    .domain([0, d3.max(hist.toArray(), d => d.count)])
    .nice();

var area = d3.area()
    .x(d => x(d.mean))
    .y0(height)
    .y1(d => y(d.count));

var blue = d3.color("steelblue"),
    darkblue = blue.darker();

// Main svg document
// Useful to place in a div with the svg in it...
// Since we're appending a group to it, we'll need to plot the parentNode later
var svg = d3.select(doc.body)
    .append("div")
        .attr("id", "viz")
    .append("svg")
        .attr("width", width+margin.left+margin.right)
        .attr("height", height+ margin.top+margin.bottom);

var viz = svg
    .append("g")
        .attr("transform",
              "translate("+margin.left+","+margin.top+")");

// Add axes groups
viz.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(0,"+height+")")
    .call(d3.axisBottom(x))

viz.append("g")
    .attr("class", "y axis")
    .call(d3.axisLeft(y));

// Add text labels...
viz.append("text")
    .attr("text-anchor", "end")
    .attr("transform", "translate("+(margin.left/2)+",0)rotate(-90)")
    .text("Count");

viz.append("text")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "central")
    .attr("transform", "translate("+(width/2)+","+(height+margin.bottom/2)+")")
        .text("Means");

// Finally, add the area element...
viz.append("path")
      .datum(hist.toArray())
      .style("stroke", darkblue)
      .style("fill", blue)
      .attr("class", "area")
      .attr("d", area);

viz.selectAll("lines")
    .data(hist.toArray())
    .enter().append("line")
        .style("stroke", darkblue)
        .attr("x1", d => x(d.mean))
        .attr("x2", d => x(d.mean))
        .attr("y1", height)
        .attr("y2", d => y(d.count));

viz.selectAll("point")
    .data(hist.toArray())
    .enter().append("circle")
        .style("fill", darkblue)
        .attr("cx", d => x(d.mean))
        .attr("cy", d => y(d.count))
        .attr("r", 2);

$$svg$$ = svg.node().parentNode.outerHTML
