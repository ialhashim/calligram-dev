#-------------------------------------------------
#
# Project created by QtCreator 2014-11-06T23:44:55
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Deformer2D
TEMPLATE = app

# Eigen
INCLUDEPATH += ../archive

SOURCES += main.cpp\
        mainwindow.cpp \
    Viewer.cpp

HEADERS  += mainwindow.h \
    Viewer.h \
    Mesh.h

FORMS    += mainwindow.ui \
    Viewer.ui
