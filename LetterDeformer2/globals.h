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

		if (segment != hole)
		{
			for (QPolygonF::iterator p = hole.begin(); p != hole.end(); p++)
			{
				std::vector< double > point;
				point.push_back(p->rx() * scaling);
				point.push_back(-p->ry() * scaling);
				holes.push_back(point);
			}
		}

        if(holes.size()) holes.resize(holes.size()-1); // remove last point, same as first
    }
};

/*
Return a RGB color value given a scalar v in the range [vmin,vmax]
In this case each color component ranges from 0 (no contribution) to
1 (fully saturated), modifications for other ranges is trivial.
The color is clipped at the end of the scales if v is outside
the range [vmin,vmax] - from StackOverflow/7706339
*/
#include <QColor>
static inline QColor qtJetColor(double v, double vmin = 0, double vmax = 1)
{
	double dv;
	if (v < vmin) v = vmin;
	if (v > vmax) v = vmax;
	dv = vmax - vmin;
	double r = 1, g = 1, b = 1;
	if (v < (vmin + 0.25 * dv)) {
		r = 0;
		g = 4 * (v - vmin) / dv;
	}
	else if (v < (vmin + 0.5 * dv)) {
		r = 0;
		b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
	}
	else if (v < (vmin + 0.75 * dv)) {
		r = 4 * (v - vmin - 0.5 * dv) / dv;
		b = 0;
	}
	else {
		g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
		b = 0;
	}
	return QColor::fromRgbF(qMax(0.0, qMin(r, 1.0)), qMax(0.0, qMin(g, 1.0)), qMax(0.0, qMin(b, 1.0)));
}
