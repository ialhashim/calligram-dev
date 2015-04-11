#pragma once
#include <vector>
#include <array>

template<typename Scalar, typename ImageClass>
std::vector< std::vector<Scalar> > make_grid(const ImageClass & img)
{
    int width = img.width(), height = img.height();
    std::vector< std::vector<Scalar> > grid(height, std::vector<Scalar>(width));
    for(int x = 0; x < width; x++)
        for(int y = 0; y < height; y++)
            grid[y][x] = (Scalar(qRed(img.pixel(x,y))) / 255.0); // change to approprate call
    return grid;
}

// 4-way recursive flooding
template<typename Scalar>
void floodfill(int x, int y, std::vector< std::vector<Scalar> > & grid, Scalar color = 0)
{
    // grid has data, and coordinates are within bound
    if(grid.empty() || grid.back().empty()) return;
    if(x < 0 || y < 0) return;
    if(x >(int)grid.at(0).size() - 1 || y >(int)grid.size() - 1) return;

    if(grid[y][x] != color){
        grid[y][x] = color;
        floodfill(x+1,y,grid,color);
        floodfill(x-1,y,grid,color);
        floodfill(x,y+1,grid,color);
        floodfill(x,y-1,grid,color);
    }
}

// 4-way stack flooding
#include <stack>
template<typename Scalar>
void floodfill_stacked(int x, int y, std::vector< std::vector<Scalar> > & grid, Scalar color = 0)
{
	if (grid.empty() || grid.back().empty()) return;

	std::stack < std::array<double, 2> > stack;
	auto makeCoord = [&](int x, int y){ std::array<double, 2> c; c[0] = x; c[1] = y;  return c; };
	stack.push(makeCoord(x,y));

	auto isFilled = [&grid,&color](int x, int y)
	{ 
		// grid has data, and coordinates are within bound
		if (x < 0 || y < 0) return true;
		if (x >(int)grid.at(0).size() - 1 || y >(int)grid.size() - 1) return true;

		return grid[y][x] == color; 
	};

	while (!stack.empty())
	{
		auto c = stack.top();
		stack.pop();
		x = c[0]; 
		y = c[1];

		if (grid[y][x] != color){
			grid[y][x] = color;
			if (!isFilled(x + 1, y)) stack.push(makeCoord(x + 1, y));
			if (!isFilled(x - 1, y)) stack.push(makeCoord(x - 1, y));
			if (!isFilled(x, y + 1)) stack.push(makeCoord(x, y + 1));
			if (!isFilled(x, y - 1)) stack.push(makeCoord(x, y - 1));
		}
	}
}

template<typename Scalar>
std::array<size_t,2> first_nonempty(std::vector< std::vector<Scalar> > & grid, Scalar empty_color = 1)
{
    std::array<size_t,2> coord = {-1, -1};
    std::array<size_t,2> sizes = {grid.size(), grid.empty() ? 0 : grid.back().size()};
    for(size_t y = 0; y < sizes[0]; y++){
        for(size_t x = 0; x < sizes[1]; x++){
            if(grid[y][x] == empty_color) continue;
            else{
                coord[0] = x;
                coord[1] = y;
                x = y = sizes[0] * sizes[1]; // break out of loops
            }
        }
    }
    return coord;
}

// Empty nighbours for coordinate (x,y) on grid
template<typename Scalar>
std::vector< std::array<size_t,2> > empty_nighbours(std::vector< std::vector<Scalar> > & grid,
                                                    const std::array<size_t,2>& cell, Scalar empty_color = 1)
{
    std::vector< std::array<size_t,2> > coords;

    for(int dy = -1; dy < 2; dy++){
        for(int dx = -1; dx < 2; dx++){
            std::array<size_t,2> c = {cell[0] + dx, cell[1] + dy};
            if(grid[c[1]][c[0]] == empty_color){
                coords.push_back(c);
            }
        }
    }

    return coords;
}

/* Adapted from https://github.com/sakri/MarchingSquaresJS */
template<typename PixelType, typename ImageType>
struct MarchingSquares
{
    enum Direction{ NONE, UP, LEFT, DOWN, RIGHT };

    ImageType image;
    double fillValue;
    int width, height;
    Direction previousStep, nextStep;

    // Change these when needed
    static double getImageValue(ImageType & image, int y, int x) {
		return image[y][x]; 
	}
    PixelType PixelTypeMaker(int x, int y){
        PixelType c;
        c[0] = x;
        c[1] = y;
        return c;
    }

    MarchingSquares( const ImageType & image, double fillValue ) : image(image), fillValue(fillValue)
    {
        width = image.front().size();
        height = image.size();
    }

    std::vector<PixelType> march()
    {
        auto coord = getFirstBorderPixelTopDown();
        if(coord[0] == -1) return std::vector<PixelType>();
        return walkPerimeter( coord );
    }

    static std::vector<PixelType> march( ImageType image, double fillValue ){
        MarchingSquares mc(image, fillValue);
        return mc.march();
    }

    PixelType getFirstBorderPixelTopDown()
    {
        PixelType pixel;
		pixel[0] = pixel[1] = -1;

        for(int y = 0; y < height; y++){
            for(int x = 0; x < width; x++){
				if (getImageValue(image, y, x) == fillValue)
				{
					pixel[0] = x;
					pixel[1] = y;
					return pixel;
				}
            }
        }

        return pixel;
    }

    std::vector<PixelType> walkPerimeter(PixelType startingPoint)
    {
        int startX = startingPoint[0];
        int startY = startingPoint[1];

        // Set up our return list
        std::vector<PixelType> pointList;

        if(startX < 0 || startY < 0) return pointList;

        // Our current x and y positions, initialized
        // to the init values passed in
        int x = startX;
        int y = startY;

        // The main while loop, continues stepping until
        // we return to our initial points
        do{
            // Evaluate our state, and set up our next direction
            step( PixelTypeMaker(x - 1, y - 1) );

            // If our current point is within our image add it to the list of points
            if (x >= 0 && x < width && y >= 0 && y < height)
                pointList.push_back( PixelTypeMaker(x - 1, y - 1) );

            switch ( nextStep )
            {
                case UP:    y--; break;
                case LEFT:  x--; break;
                case DOWN:  y++; break;
                case RIGHT: x++; break;
            }

        } while (x != startX || y != startY);

        pointList.push_back( PixelTypeMaker(x - 1, y - 1) );

        return pointList;
    }

    inline bool emptyPixel( int x, int y )
    {
        if(x < 0 || x > width - 1 || y < 0 || y > height - 1) return true;
        return getImageValue(image,y,x) != fillValue;
    }

    void step( PixelType pixel )
    {
        bool upLeft = emptyPixel( pixel[0], pixel[1] );
        bool upRight = emptyPixel( pixel[0]+1, pixel[1] );
        bool downLeft = emptyPixel( pixel[0], pixel[1]+1 );
        bool downRight = emptyPixel( pixel[0]+1, pixel[1]+1 );

        // Store our previous step
        previousStep = nextStep;

        // Determine which state we are in
        int state = 0;
        if (upLeft)		state |= 1;
        if (upRight)	state |= 2;
        if (downLeft)	state |= 4;
        if (downRight)	state |= 8;

        // So we can use a switch statement to determine our
        // next direction based on
        switch (state){
        case 1: nextStep = UP; break;
        case 2: nextStep = RIGHT; break;
        case 3: nextStep = RIGHT; break;
        case 4: nextStep = LEFT; break;
        case 5: nextStep = UP; break;
        case 6:
            if (previousStep == UP){
                nextStep = LEFT;
            }else{
                nextStep = RIGHT;
            }
            break;
        case 7: nextStep = RIGHT; break;
        case 8: nextStep = DOWN; break;
        case 9:
            if (previousStep == RIGHT){
                nextStep = UP;
            }else{
                nextStep = DOWN;
            }
            break;
        case 10: nextStep = DOWN; break;
        case 11: nextStep = DOWN; break;
        case 12: nextStep = LEFT; break;
        case 13: nextStep = UP; break;
        case 14: nextStep = LEFT; break;
        default:
            nextStep = NONE;//this should never happen
            break;
        }
    }
};

// Qt specific utility functions:
QPolygonF resamplePolygon(QPolygonF points, int count = 100){
    QPainterPath path;
    path.addPolygon(points);
    auto pathLen = path.length();
    auto stepSize = pathLen / count;
    QPolygonF newPoints;
    for(int i = 0; i < count; i++)
        newPoints << path.pointAtPercent(path.percentAtLength( stepSize * i ));
    newPoints << path.pointAtPercent(0);
    return newPoints;
}

QPolygonF smoothPolygon( QPolygonF points, int iterations = 1 ){
    for(int i = 0; i < iterations; i++)
    {
        QPolygonF newPoints;
        points.removeLast();

        for(int p = 0; p < points.size(); p++){
            int s = p-1, t = p+1;
            if(s < 0) s = points.size() - 1;
            if(t > points.size()-1) t = 0;
            QPointF prev = points[s];
            QPointF next = points[t];
            newPoints << (prev + next) * 0.5;
        }

        newPoints << newPoints.front();
        points = newPoints;
    }
    return points;
}
