QT -= gui

TARGET = libShapeOp
TEMPLATE = lib
CONFIG += staticlib

INCLUDEPATH += .

SOURCES += API.cpp
HEADERS += API.h \
    Common.h \
    Constraint.h \
    Force.h \
    LSSolver.h \
    Solver.h \
    Types.h

# Build options
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}
DESTDIR = $$PWD/$$CFG/lib

win32:DEFINES += SHAPEOP_MSVC SHAPEOP_OPENMP SHAPEOP_EXPORT SHAPEOP_HEADER_ONLY
