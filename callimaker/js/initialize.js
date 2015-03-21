// Background color
$("#content").css("background-color", randomBackgroundColor());

// Menu button animation:
$("#menu-toggle-button").click(function() {
    if(menuState.state === "hidden")
    {
        TweenLite.to("#menu-toggle-button", 0.15, {backgroundColor:"hsl(0,0%,77%)"} );
        TweenLite.to("#menu-button-icon", 0.15, {rotation:45});

        TweenLite.to("#menu-wrapper", 0.25, {top:"0px", ease: Back.easeOut.config(1), delay: 0.15});
        menuState.state = "shown";
    }
    else
    {
        TweenLite.to("#menu-toggle-button", 0.15, {backgroundColor:"hsl(0,83%,77%)"} );
        TweenLite.to("#menu-button-icon", 0.15, {rotation:0});

        TweenLite.to("#menu-wrapper", 0.3, {top:"-150px", delay: 0.15});
        menuState.state = "hidden";
    }
});

// Load example shape:
selectedShape.name = "kangaroo";
selectedShape.geometry = exampleShape;

$("#selected-word").text(selectedShape.name);
$("#active-shape").attr("d", selectedShape.geometry );


// Initialize drawing area:
userCanvas = document.getElementById('user-canvas');
userCanvas.context = userCanvas.getContext("2d");

// Simple drawing code:
var context = userCanvas.context;
var clickX = new Array;
var clickY = new Array;
var clickDrag = new Array;
var paint;

function addClick(x, y, dragging)
{
  clickX.push(x);
  clickY.push(y);
  clickDrag.push(dragging);
}

function redraw(){
  context.clearRect(0, 0, context.canvas.width, context.canvas.height);
  context.strokeStyle = "white";
  context.lineJoin = "round";
  context.lineWidth = 6;

  for(var i=0; i < clickX.length; i++) {
    context.beginPath();
    if(clickDrag[i] && i){
      context.moveTo(clickX[i-1], clickY[i-1]);
     }else{
       context.moveTo(clickX[i]-1, clickY[i]);
     }
     context.lineTo(clickX[i], clickY[i]);
     context.closePath();
     context.stroke();
  }
}

$('#user-canvas').mousedown(function(e){
  var mouseX = e.pageX - this.offsetLeft;
  var mouseY = e.pageY - this.offsetTop;

  paint = true;
  addClick(e.pageX - this.offsetLeft, e.pageY - this.offsetTop);
  redraw();
});

$('#user-canvas').mousedown(function(e){
  var mouseX = e.pageX - this.offsetLeft;
  var mouseY = e.pageY - this.offsetTop;

  paint = true;
  addClick(e.pageX - this.offsetLeft, e.pageY - this.offsetTop);
  redraw();
});

$('#user-canvas').mousemove(function(e){
  if(paint){
    addClick(e.pageX - this.offsetLeft, e.pageY - this.offsetTop, true);
    redraw();
  }
});

$('#user-canvas').mouseup(function(e){
  paint = false;
  if (e.which === 3)
  {
      clickX = new Array;
      clickY = new Array;
      clickDrag = new Array;
  }
});

$('#user-canvas').mouseleave(function(e){
  paint = false;
});

