QT       += core gui opengl
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET       = tool
TEMPLATE     = app
INCLUDEPATH += . ..

SOURCES  += main.cpp  \
            mainwindow.cpp \
            Polygon.cpp \
            Viewer.cpp

HEADERS  += mainwindow.h \
            Polygon.h \
            Viewer.h

FORMS    += mainwindow.ui \
            Viewer.ui

### Libraries:

# Triangle library
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}
LIBS += -L$$PWD/../triangle/$$CFG/lib -ltriangle
INCLUDEPATH += ../triangle
DEFINES *= NO_TIMER ANSI_DECLARATORS TRILIBRARY

# Voronoi library
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}
LIBS += -L$$PWD/../voronoi/$$CFG/lib -lvoronoi
INCLUDEPATH += ../voronoi
