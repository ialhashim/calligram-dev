#include "MainWindow.h"
#include "ui_MainWindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    Q_INIT_RESOURCE(MainWindow);

    ui->setupUi(this);

    ui->webView->load(QUrl("qrc:/view.html"));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::keyPressEvent(QKeyEvent *event)
{
    switch(event->key())
    {
        case Qt::Key_Escape: close();break;
        default: QMainWindow::keyPressEvent(event);
    }
}
