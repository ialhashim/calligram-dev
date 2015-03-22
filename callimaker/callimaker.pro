#-------------------------------------------------
#
# Project created by QtCreator 2015-03-19T23:37:16
#
#-------------------------------------------------

QT       += core gui webkitwidgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = callimaker
TEMPLATE = app


SOURCES += main.cpp\
        MainWindow.cpp \
    Control.cpp \
    WebPage.cpp

HEADERS  += MainWindow.h \
    Control.h \
    WebPage.h

FORMS    += MainWindow.ui

RC_FILE = MainWindow.rc

RESOURCES += MainWindow.qrc
