> Repository: https://github.com/alangrafu/icicle-chart-d3 <br>
> Author: Alvaro Graves <br>
> License: Apache-2.0 license, see [LICENSE_icicle-chart](./LICENSE_icicle-chart) file in this folder,

# Icicle Chart
A reusable icicle chart implementation in D3.js. Styleable, configurable and transition-capable.

[Icicle Example](http://bl.ocks.org/tpreusse/2bc99d74a461b8c0acb1)

## Usage

### Install
`bower install git@github.com:alangrafu/icicle-chart-d3.git --save`

### Data structure
```
var data = [
  {
    className: 'germany', // optional, can be used for styling
    axes: [
      {axis: "strength", value: 13, yOffset: 10},
      {axis: "intelligence", value: 6},
      {axis: "charisma", value: 5},  
      {axis: "dexterity", value: 9},  
      {axis: "luck", value: 2, xOffset: -20}
    ]
  },
  {
    className: 'argentina',
    axes: [
      {axis: "strength", value: 6},
      {axis: "intelligence", value: 7},
      {axis: "charisma", value: 10},  
      {axis: "dexterity", value: 13},  
      {axis: "luck", value: 9}
    ]
  }
];
```

`xOffset` and `yOffset` are optional values that allows a developer to change the position of a specific label. It is important to add them in **the first** group of axes.

### Simple single chart drawing
```html
<div class="chart-container"></div>
<script>
IcicleChart.draw(".chart-container", data);
</script>
```

### D3.js reusable chart API
```javascript
var chart = IcicleChart.chart();
var svg = d3.select('body').append('svg')
  .attr('width', 600)
  .attr('height', 800);

// draw one
svg.append('g').classed('focus', 1).datum(data).call(chart);

// draw many icicles
var game = svg.selectAll('g.game').data(
  [
    data,
    data,
    data,
    data
  ]
);
game.enter().append('g').classed('game', 1);
game
  .attr('transform', function(d, i) { return 'translate(150,600)'; })
  .call(chart);
```

### Style with CSS
```css
.icicle-chart .area {
  fill-opacity: 0.7;
}
.icicle-chart.focus .area {
  fill-opacity: 0.3;
}
.icicle-chart.focus .area.focused {
  fill-opacity: 0.9;
}
.area.germany, .germany .circle {
  fill: #FFD700;
  stroke: none;
}
.area.argentina, .argentina .circle {
  fill: #ADD8E6;
  stroke: none;
}
```

### Configure
```javascript
// retrieve config
chart.config();
// all options with default values
chart.config({
  containerClass: 'icicle-chart', // target with css, the default stylesheet targets .icicle-chart
  w: 600,
  h: 600,
  factor: 0.95,
  factorLegend: 1,
  levels: 3,
  maxValue: 0,
  minValue: 0,
  radians: 2 * Math.PI,
  color: d3.scale.category10(), // pass a noop (function() {}) to decide color via css
  axisLine: true,
  axisText: true,
  circles: true,
  radius: 5,
  open: false,  // whether or not the last axis value should connect back to the first axis value
                // if true, consider modifying the chart opacity (see "Style with CSS" section above)
  axisJoin: function(d, i) {
    return d.className || i;
  },
  tooltipFormatValue: function(d) {
    return d;
  },
  tooltipFormatClass: function(d) {
    return d;
  },
  transitionDuration: 300
});
```

## Example
### CSV2Icicle

Display a csv file as a icicle chart at [http://alangrafu.github.io/icicle-chart-d3/csv2icicle.html](http://alangrafu.github.io/icicle-chart-d3/csv2icicle.html).


[![Example](https://rawgit.com/tpreusse/icicle-chart-d3/master/example/demo.svg)](http://bl.ocks.org/tpreusse/2bc99d74a461b8c0acb1)

