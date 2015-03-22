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

// Drawing code:
$( document ).ready(function() {
    console.log( "Page ready." + Control.callableEverywhere() );

    var dp = new DrawPanel("user-canvas");
});
