#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "Control.h"
#include "WebPage.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    Q_INIT_RESOURCE(MainWindow);

    ui->setupUi(this);

    auto webpage = new WebPage(this);
    ui->webView->setPage(webpage);
    webpage->mainFrame()->load(QUrl("qrc:/view.html"));

    // Connect control
    auto frame = ui->webView->page()->mainFrame();
    auto c = new Control(frame);
    frame->addToJavaScriptWindowObject("Control", c);
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
