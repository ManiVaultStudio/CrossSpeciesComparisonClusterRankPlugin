// ManiVault invokes this function to set the plot data,
// when emitting qt_js_setDataInJS from the communication object
// The connection is established in qwebchannel.tools.js


window.onresize = doALoadOfStuff;
var data;
var svg;
var x;
var y;
var partition;
var rect;
var margin = { top: 10, right: 10, bottom: 10, left: 10 };
var width;
var height;
var selectedClusterContainer = [];
function doALoadOfStuff() {
    if (data != "") {
        VisCreate();
    }


}

function removeRectHighlight(d) {
    if (data != "") {
        d3.selectAll('rect')
            .style('stroke', '#000')
            .style('stroke-width', 0.2);
        selectedClusterContainer = [];

    }
}


function drawChart(d) {

    if (d != "") {
        //log("\nData received: \n" + d)
        var parsedData = JSON.parse(d);
        // If parsedData is an array, take the first element
        data = Array.isArray(parsedData) ? parsedData[0] : parsedData;
        VisCreate();

        /*data = {
            name: "flare",
            color: "pink",
            children: [
                {
                    name: "analytics",
                    color: "blue",
                    children: [
                        {
                            name: "cluster",
                            color: "green",
                            children: [
                                { name: "AgglomerativeCluster", color: "red", value: 3938 },
                                { name: "CommunityStructure", color: "orange", value: 3812 },
                                { name: "HierarchicalCluster", color: "yellow", value: 6714 },
                                { name: "MergeEdge", color: "purple", value: 743 },
                            ],
                        },
                        {
                            name: "graph",
                            color: "cyan",
                            children: [
                                { name: "BetweennessCentrality", color: "magenta", value: 3534 },
                                { name: "LinkDistance", color: "lime", value: 5731 },
                                { name: "MaxFlowMinCut", color: "teal", value: 7840 },
                                { name: "ShortestPaths", color: "maroon", value: 5914 },
                                { name: "SpanningTree", color: "navy", value: 3416 },
                            ],
                        },
                        {
                            name: "optimization",
                            color: "olive",
                            children: [
                                { name: "AspectRatioBanker", color: "silver", value: 7074 },
                            ],
                        },
                    ],
                },
            ],
        }; */
        VisCreate();
    }
}

