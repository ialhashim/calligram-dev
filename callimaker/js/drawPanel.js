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

            var majorStrokeLayer = new Layer();

            majorStroke = new Path();
            majorStroke.selected = false;
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
                majorStroke.strokeColor = 'rgba(255,255,255,0.75)';
                majorStroke.strokeWidth = 10;
                majorStroke.shadowColor =  'rgba(0,0,0,0.2)';
                majorStroke.shadowBlur = 8;
                majorStroke.shadowOffset = new Point(2, 2);
            }
            else
                return;

            // Convert to SVG text:
            //var majorStrokeSVG = majorStroke.exportSVG(true);

            // Draw letters at stroke:
            var lettersLayer = new Layer();

            var nc = selectedShape.name.length;
            var step = majorStroke.length / nc;
            for(var i = 0; i < nc; i++)
            {
                var c = selectedShape.name[i];
                var t = i/nc;
                var alpha = (t * 0.4) + 0.5;
                var l = i * step;
                var position = majorStroke.getPointAt(l);

                var text = new PointText({
                    point: position,
                    content: c,
                    fillColor: 'rgba(0,0,0,'+ alpha +')',
                    fontFamily: 'VAGRound',
                    fontWeight: 'bold',
                    fontSize: 120,
                    shadowColor :  'rgba(0,0,0,0.8)',
                    shadowBlur : 8,
                    shadowOffset : new Point(2, 2)
                });

                text.translate(text.globalToLocal(text.bounds.bottomLeft).add(new Point(step*0.2,0)));
            }

            // Send to back-end:
            //console.log(majorStrokeSVG);
            //var data = {command: "majorStroke", data: majorStrokeSVG};
            //Control.receiveData(data);

            // Change to letters tool:
            var lettersStrokesLayer = new Layer();
            lettersStrokeTool.activate();

            // Create default letter strokes:
            for(var i = 1; i < nc; i++)
            {
                var c = selectedShape.name[i];
                var t = i/nc;
                var alpha = (t * 0.4) + 0.5;
                var l = i * step;
                var position = majorStroke.getPointAt(l);
                var normal = majorStroke.getNormalAt(l);
                normal.length = 80;

                var from = position.add(normal);
                var to = position.add(normal.negate());

                var letterStroke = new Path.Line(from, to);
                lsColor = 'hsla('+(t * 128)+',100%,70%,1.0)';
                letterStroke.strokeColor = lsColor;
                letterStroke.strokeCap = 'round';
                letterStroke.strokeJoin = 'round';
                letterStroke.strokeWidth = 10;
                letterStroke.shadowColor =  'rgba(0,0,0,0.2)';
                letterStroke.shadowBlur = 8;
                letterStroke.shadowOffset = new Point(2, 2);
            }
        }
    }

    // 2) Letter Strokes:
    if(lettersStrokeTool)
    {
        lettersStrokeTool.onMouseDown = function(event)
        {
            var closestLS, closestDist = 1e10;

            var existingStrokes = project.activeLayer.children;
            for(var i = 0; i < existingStrokes.length; i++)
            {
                var stroke = existingStrokes[i];
                var dist = stroke.getNearestPoint(event.point).subtract(event.point).length;

                if(dist < closestDist){
                    closestLS = existingStrokes[i];
                    closestDist = dist;
                }
            }

            var lsColor = closestLS.strokeColor;
            closestLS.remove();

            letterStroke = new Path();
            letterStroke.strokeColor = 'rgba(255,255,0,0.50)';
            letterStroke.strokeColor = lsColor;
            letterStroke.strokeCap = 'round';
            letterStroke.strokeJoin = 'round';
            letterStroke.strokeWidth = 10;
            letterStroke.shadowColor =  'rgba(0,0,0,0.2)';
            letterStroke.shadowBlur = 8;
            letterStroke.shadowOffset = new Point(2, 2);
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
