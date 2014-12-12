QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = IsometricWillmoreFlow1D
TEMPLATE = app

INCLUDEPATH += .

SOURCES += main.cpp \
    IsometricWillmoreFlow1D.cpp \
    Mesh.cpp \
    Vertex.cpp

HEADERS += \
    IsometricWillmoreFlow1D.h \
    Mesh.h \
    Vertex.h \
    HalfEdge.h
