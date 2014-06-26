#include "Viewer.h"
#include "ui_Viewer.h"

Viewer::Viewer(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Viewer)
{
    ui->setupUi(this);
}

Viewer::~Viewer()
{
    delete ui;
}

void Viewer::paintEvent(QPaintEvent *)
{

}
