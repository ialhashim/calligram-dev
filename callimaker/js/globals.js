var background_colors = ["#84A5DD","#64d6e2","#FDAED4","#785ebb","#a09de5",
                         "#A0CADB","#79BBB5","#F6D860","#9988CD","#EDB948",
                         "#FDACB4","#dbbe39","#DFBC94","#3CCAD1","#F68F6F",
                         "#9787EA","#ccc5e3","#FD9372","#9B7AD5","#de89ac",
                         "#98BFF6","#FEC54F","#9784ED","#F4C3C5","#63C5AB",
                         "#F4696B","#E794AE","#9B7FE6","#EF9F64","#87C4A3"];

function randomBackgroundColor(){
	return background_colors[Math.floor(Math.random()*background_colors.length)];
}

var menuState = {state: "hidden"};

var userStroke = {};
