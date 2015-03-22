paper.install(window);

function DrawPanel( canvasID )
{
    var canvas = $("#" + canvasID)[0];
    var ctx = canvas.getContext('2d');

    paper.setup(canvas);

    /// Drawing tools:

    var majorStrokeTool = new Tool();
    var lettersStrokeTool = new Tool();

    var majorStroke;
    var letterStroke;

    // 1) Major Stroke:
    if(majorStrokeTool)
    {
        // Define a mousedown and mousedrag handler
        majorStrokeTool.onMouseDown = function(event)
        {
            if(majorStroke) majorStroke.remove();

            majorStroke = new Path();
            majorStroke.selected = true;
            majorStroke.strokeColor = 'white';
            majorStroke.strokeCap = 'round';
            majorStroke.strokeJoin = 'round';
            majorStroke.strokeWidth = 20;
            majorStroke.add(event.point);
        }

        majorStrokeTool.onMouseDrag = function(event) {
            majorStroke.add(event.point);
        }

        majorStrokeTool.onMouseUp = function(event)
        {
            // Post process:
            if(majorStroke.length !== 0)
            {
                majorStroke.setSegments( pathUtils.laplacianSmoothing(majorStroke, 10) );
                majorStroke.simplify();
            }
            else
                return;

            // Convert to SVG text:
            var majorStrokeSVG = majorStroke.exportSVG(true);

            // Draw letters at stroke:
            var nc = selectedShape.name.length;
            var step = majorStroke.length / nc;
            for(var i = 0; i < nc; i++)
            {
                var c = selectedShape.name[i];
                var text = new PointText({
                    point: majorStroke.getPointAt(i * step),
                    content: c,
                    fillColor: 'black',
                    fontFamily: 'VAGRound',
                    fontWeight: 'bold',
                    fontSize: 120
                });
            }

            // Send to back-end:
            console.log(majorStrokeSVG);
            var data = {command: "majorStroke", data: majorStrokeSVG};
            Control.receiveData(data);

            // Change to letters tool:
            lettersStrokeTool.activate();
        }
    }

    // 2) Letter Strokes:
    if(lettersStrokeTool)
    {
        lettersStrokeTool.onMouseDown = function(event)
        {
            letterStroke = new Path();
            letterStroke.strokeColor = 'rgba(255,255,0,0.50)';
            letterStroke.strokeCap = 'round';
            letterStroke.strokeJoin = 'round';
            letterStroke.strokeWidth = 10;
            letterStroke.add(event.point);
        }

        lettersStrokeTool.onMouseDrag = function(event) {
            letterStroke.add(event.point);
        }

        lettersStrokeTool.onMouseUp = function(event){
            if(letterStroke.length !== 0)
            {
                letterStroke.setSegments( pathUtils.laplacianSmoothing(letterStroke, 10) );
                letterStroke.simplify();
            }
            else
                return;
        }
    }

    console.log( "Drawing panel ready." );
}

var Utils = function(){};
Utils.prototype.laplacianSmoothing = function (orig_path, iterations){
    if(!orig_path._segments.length) return orig_path._segments;
    if(orig_path.length === 0) return orig_path._segments;

    // Resample points of path
    var points = this.points = [];
    var numSteps = 50;
    var step = orig_path.length / numSteps;
    for(var i = 0; i < numSteps; i++)
    {
        offset = i * step;
        points.push(orig_path.getPointAt(offset));
    }

    // Apply laplacian smoothing on points
    for (var itr = 0; itr < iterations; itr++)
    {
        var smoothPoints = [];
        for (var i = 0; i < points.length; i++)
        {
            if(i === 0 || i === points.length - 1)
                smoothPoints.push(points[i]);
            else
                smoothPoints.push(points[i-1].add(points[i+1]).divide(2).clone());
        }

        // Clone results
        points = [];
        for (var i = 0; i < smoothPoints.length; i++)
            points.push(smoothPoints[i].clone());
    }

    // Create segments from path
    var newsegments = [];
    for(var i = 0; i < points.length - 1; i++)
        newsegments.push(new Segment(points[i]), new Segment(points[i+1]));
    return newsegments;
}
var pathUtils = new Utils;
