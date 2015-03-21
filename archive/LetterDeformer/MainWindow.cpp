#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "Viewer.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->mainlayout->addWidget(new Viewer());
}

MainWindow::~MainWindow()
{
    delete ui;
}
