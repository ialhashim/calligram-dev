#pragma once

#include <QGLWidget>
#include <QVector>

namespace surface_mesh{class Mesh;}

namespace Ui {
class Viewer;
}

class MyQImage : public QImage{public: int x,y; double opacity; MyQImage(const QImage & img = QImage(),
                               int x = 0, int y = 0, double opacity = 1.0) : QImage(img), x(x), y(y), opacity(opacity) {}};

class Viewer : public QGLWidget
{
    Q_OBJECT

public:
    explicit Viewer(QWidget *parent = 0);
    ~Viewer();

    void paintEvent(QPaintEvent *);
    void keyPressEvent(QKeyEvent *event);
    void mousePressEvent(QMouseEvent *event);

    QVector<MyQImage> images;

	QVector<surface_mesh::Mesh*> meshes;	   
	
	// Cage
	QPolygonF points;
	QVector<double> control_edge_length;
	QVector<QPointF> control_edge_normal;
	bool isCoordinatesReady;

	QPolygonF chull;

private:
    Ui::Viewer *ui;

public slots:
    void demo();

	// code from: https://github.com/Juyong/OpenFlipper/Plugin-LBC
	void computeGreenCoordinates();
	void applyDeformGreen();
	void applyDeformMVC();
};

// Utility:
#pragma once
#include <deque>

// Melkman's algorithm for finding convex hull of 2D points
// From CGAL
template <class Point, class InputIterator, class OutputIterator>
static inline OutputIterator convexhull2d(InputIterator first, InputIterator last, OutputIterator result)
{
    // Subroutine
    auto left_turn = [&](const Point& a, const Point& b, const Point& c){
        return ((b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1]) > 0);
    };
    auto equal_points = [&](const Point& a, const Point& b){
      return a[0] == b[0] && a[1] == b[1];
    };

    std::deque<Point> Q;
    if (first == last) return result;								// 0 elements
    Point p = *first;
    if (++first == last) { *result = p; ++result; return result; }	// 1 element
    Point q = *first;
    if (++first == last) {											// 2 elements
        *result = p; ++result;
        if (!equal_points(p, q)){ *result = q; ++result; }
        return result;
    }
    Q.push_back(p);

    Point r;
    while (first != last){
        r = *first;
        // visited input sequence =  p,..., q, r
        if (left_turn(p, q, r)) { Q.push_back(q);  break; }
        if (left_turn(q, p, r)) { Q.push_front(q); break; }
        q = r;
        ++first;
    }

    Point current = q;
    if (first != last)           // non-collinear point r exists
    {
        current = r;
        // current, Q.front(), ..., Q.back()
        // ccw convex hull of visited points
        Point s;
        while (++first != last)
        {
            r = *first;
            if (left_turn(current, r, Q.front()) || left_turn(Q.back(), r, current))
            {
                s = current;
                while (!Q.empty() && !left_turn(r, s, Q.front())) {s = Q.front(); Q.pop_front();}
                Q.push_front(s);
                s = current;
                while (!Q.empty() && !left_turn(s, r, Q.back())) {s = Q.back(); Q.pop_back();}
                Q.push_back(s);
                current = r;
            }
        }
    }

    Q.push_back(current);
    std::copy(Q.begin(), Q.end(), result);
    return result;
}
