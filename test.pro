#-------------------------------------------------
#
# Project created by QtCreator 2018-03-14T16:38:37
#
#-------------------------------------------------

QT       += core gui
QT       += opengl widgets
QT += opengl
LIBS = -lGLU \

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = test
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    Window.cpp \
    function.cpp \
    intersection_func.cpp \
    predicates.c \
    ttotr.c

HEADERS  += mainwindow.h \
    Window.h \
    geom_func.h \
    datastruct.h

FORMS    += mainwindow.ui

DISTFILES += \
    surfacemesh.vtk \
    spacemesh.vtk
