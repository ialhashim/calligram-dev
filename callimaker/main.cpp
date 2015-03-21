#include "MainWindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    //w.setWindowFlags( Qt::Window | Qt::WindowTitleHint );
    w.show();

    return a.exec();
}
