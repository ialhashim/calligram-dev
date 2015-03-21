TEMPLATE = app

INCLUDEPATH += .

HEADERS += \
    allocMultidim.h \
    GLWindow.h \
    Grid.h \
    pngLoad.h \
    Vector2.h \
    matrix.hh \
    MatrixStripped.hh

SOURCES += \
    GLWindow.cpp \
    Grid.cpp \
    main.cpp \
    pngLoad.cpp \
    allocMultidim.cpp

# Copy freeglut.dll
win32{
    LIBS += $$PWD/GL/freeglut.lib

    EXTRA_BINFILES += $$PWD/GL/freeglut.dll

    EXTRA_BINFILES_WIN = $${EXTRA_BINFILES}
    EXTRA_BINFILES_WIN ~= s,/,\\,g
    DESTDIR_WIN = $${EXECUTABLEPATH}
    DESTDIR_WIN ~= s,/,\\,g
    for(FILE,EXTRA_BINFILES_WIN){
        # message("Will copy file" $$FILE "to" $$DESTDIR_WIN)
        QMAKE_POST_LINK += $$quote(cmd /c copy /y $${FILE} $${DESTDIR_WIN}$$escape_expand(\n\t))
    }
}
