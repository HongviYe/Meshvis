#include "mainwindow.h"
#include <QApplication>
#include <QLabel>
#include <QPushButton>
#include <QDebug>
#include <QtOpenGL>
#include <QWidget>
#include <qmessagebox.h>
#include <Window.h>
#include <GL/glu.h>
#include <datastruct.h>
void myDisplay(void)

{

     glClear(GL_COLOR_BUFFER_BIT);

     glRectf(-0.5f, -0.5f, 0.5f, 0.5f);

     glFlush();

}

int main(int argc, char *argv[])
{
    bool fs = false;
    char *fname1="/home/yhf/test/spacemesh.vtk";
    char *fname2="/home/yhf/test/surfacemesh.vtk";
    mesh F6(fname1,fname2);
    QApplication app(argc, argv);
    /*
    switch( QMessageBox::information( 0,
          "Start FullScreen?",
          "Would You Like To Run In Fullscreen Mode?",
          QMessageBox::Yes,
          QMessageBox::No | QMessageBox::Default ) )
      {
      case QMessageBox::Yes:
        fs = true;
        break;
      case QMessageBox::No:
        fs = false;
        break;
      }
    //*/
    QSurfaceFormat format;
        format.setSamples(4);
    Window window;
    window.flight=&F6;
    window.setFormat(format);
    window.resizeGL(600,600);
   // std::cout<<F6.BL.size()<<std::endl;
    window.show();
    //MainWindow a;
    //a.setWindowTitle("woshi nibaba");
    //a.show();


        return app.exec();


}
