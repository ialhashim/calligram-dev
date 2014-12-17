#pragma once

#include <QPainterPath>

struct LetterGeometry{
    std::vector< std::vector< double > > points, holes;
    LetterGeometry(QChar c = 'n', int fontsize = 40, double scaling = 0.3){
        QPainterPath text;
        text.addText(0,0,QFont("VAG Rounded", fontsize), c);
        QList<QPolygonF> poly = text.toSubpathPolygons();

        auto segment = poly.front();
        auto hole = poly.back();

        for (QPolygonF::iterator p = segment.begin(); p != segment.end(); p++)
        {
            std::vector< double > point;
            point.push_back(p->rx() * scaling);
            point.push_back(-p->ry() * scaling);
            points.push_back(point);
        }
        points.resize(points.size()-1); // remove last point, same as first

        for (QPolygonF::iterator p = hole.begin(); p != hole.end(); p++)
        {
            std::vector< double > point;
            point.push_back(p->rx() * scaling);
            point.push_back(-p->ry() * scaling);
            holes.push_back(point);
        }
        holes.resize(holes.size()-1); // remove last point, same as first
    }
};
