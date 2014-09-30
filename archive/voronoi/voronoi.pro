TARGET = voronoi
TEMPLATE = lib
CONFIG += staticlib

SOURCES += voronoi.cpp
HEADERS += voronoi.h

# Build flag
CONFIG(debug, debug|release) {CFG = debug} else {CFG = release}

# LIB output folder
DESTDIR = $$PWD/$$CFG/lib
