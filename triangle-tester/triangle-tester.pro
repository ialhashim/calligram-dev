QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET   = triangle-tester
TEMPLATE = app

SOURCES  += main.cpp mainwindow.cpp
HEADERS  += mainwindow.h

FORMS    += mainwindow.ui

# Libraries
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}
LIBS += -L$$PWD/../triangle/$$CFG/lib -ltriangle
INCLUDEPATH += ../triangle
DEFINES *= NO_TIMER ANSI_DECLARATORS TRILIBRARY
