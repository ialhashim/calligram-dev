QT       += core gui opengl
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = HandleHoles
TEMPLATE = app

SOURCES  += main.cpp mainwindow.cpp \
            Viewer.cpp
HEADERS  += mainwindow.h \
            Viewer.h
FORMS    += mainwindow.ui \
            Viewer.ui
