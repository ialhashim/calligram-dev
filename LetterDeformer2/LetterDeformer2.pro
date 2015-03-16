QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = LetterDeformer
TEMPLATE = app

SOURCES     +=  main.cpp MainWindow.cpp \
                Viewer.cpp
HEADERS     +=  MainWindow.h \
                Viewer.h \
                Mesh.h \
    globals.h
FORMS       +=  MainWindow.ui


# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

## Libraries:
# ShapeOp
LIBS += -L$$PWD/libShapeOp/$$CFG/lib -llibShapeOp
INCLUDEPATH += ./libShapeOp
DEFINES += SHAPEOP_EXPORT

# Triangle
LIBS += -L$$PWD/libTriangle/$$CFG/lib -ltriangle
INCLUDEPATH += ./libTriangle
DEFINES += NO_TIMER ANSI_DECLARATORS TRILIBRARY _CRT_SECURE_NO_WARNINGS

win32-msvc*{
	QMAKE_CXXFLAGS *= /openmp
}
