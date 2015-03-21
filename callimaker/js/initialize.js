// Background color
$("#content").css("background-color", randomBackgroundColor());

// Menu button
$("#menu-toggle-button").click(function() {

    if(menuState.state === "hidden")
    {
        TweenLite.to("#menu-wrapper", 0.25, {top:"0px", ease: Back.easeOut.config(1)});
        menuState.state = "shown";
    }
    else
    {
        TweenLite.to("#menu-wrapper", 0.3, {top:"-150px"});
        menuState.state = "hidden";
    }
});
