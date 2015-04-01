#-------------------------------------------------
#
# Project created by QtCreator 2015-03-31T21:49:53
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = SimpleFitting
TEMPLATE = app


SOURCES +=  main.cpp\
            mainwindow.cpp \
            Viewer.cpp

HEADERS  += mainwindow.h \
            Viewer.h

FORMS    += mainwindow.ui \
            Viewer.ui


# Build flag
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

# Libraries
INCLUDEPATH += $$PWD/../libfastmarching
LIBS += $$PWD/../libfastmarching/$$CFG/lib/libfastmarching.lib
