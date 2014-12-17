TARGET = triangle
TEMPLATE = lib
CONFIG += staticlib

SOURCES += triangle.cpp
HEADERS += triangle.h trianglelib.h

# Build flag
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

# LIB output folder
DESTDIR = $$PWD/$$CFG/lib

DEFINES += NO_TIMER ANSI_DECLARATORS TRILIBRARY _CRT_SECURE_NO_WARNINGS
