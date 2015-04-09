#pragma once

int globalStrokeSize = 25;

template<typename QPolygonType>
QPolygonType resamplePolygon(QPolygonType points, int count = 100){
    QPainterPath path;
    path.addPolygon(points);
    auto pathLen = path.length();
    auto stepSize = pathLen / count;
    QPolygonType newPoints;
    for(int i = 0; i < count; i++)
    {
        QPolygon::value_type p = path.pointAtPercent(path.percentAtLength( stepSize * i )).toPoint();
        newPoints << p;
    }
    return newPoints;
}

template<typename QPolygonType>
QPolygonType smoothPolygon( QPolygonType points, int iterations = 1 ){
    for(int i = 0; i < iterations; i++)
    {
        QPolygonType newPoints;

        for(int p = 0; p < points.size(); p++)
        {
            int s = p-1, t = p+1;
            if(s < 0) s = 0;
            if(t > points.size()-1) t = points.size()-1;
            auto prev = points[s];
            auto next = points[t];
            newPoints << (prev + next) * 0.5;
        }

        points = newPoints;
    }
    return points;
}

MyQImage visualizeStroke(QImage shape, QPolygon path)
{
    QImage strokeImage(shape.width(), shape.height(), QImage::Format_ARGB32_Premultiplied);
    strokeImage.fill(Qt::transparent);
    QPainter painter(&strokeImage);
    painter.setRenderHint(QPainter::Antialiasing);

    // Stroke color
    QColor penColor;
    penColor.setHslF(double(rand()) / RAND_MAX * 0.25, 1, 0.5);
    painter.setPen(QPen(penColor, globalStrokeSize, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

    // Draw stroke
    painter.drawPolyline(path);

    return strokeImage;
}
