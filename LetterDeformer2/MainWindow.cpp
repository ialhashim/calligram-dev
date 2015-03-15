#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "Viewer.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
	auto v = new Viewer();
    ui->mainlayout->addWidget(v);

	v->setFocus(Qt::NoFocusReason);
}

MainWindow::~MainWindow()
{
    delete ui;
}
