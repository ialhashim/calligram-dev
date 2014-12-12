#include <QApplication>

#include "IsometricWillmoreFlow1D.h"
#include <QSharedPointer>
#include "Viewer.h"

int main(int argc, char *argv[])
{
	if (argc < 2) return -1;

    auto flow = QSharedPointer<IsometricWillmoreFlow1D>(new IsometricWillmoreFlow1D(new Mesh(argv[1])));

    QApplication app(argc, argv);

    auto viewer = new Viewer(flow);
    viewer->show();

    return app.exec();
}
