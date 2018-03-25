#ifndef WINDOW_H
#define WINDOW_H
#include <QOpenGLWindow>
#include <QWidget>
#include <datastruct.h>
#include <QGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLWidget>
#include <QtOpenGL>
#include <cmath>
#include <QOpenGLBuffer>
#include <QOpenGLContext>
#include <QOpenGLShader>
#include "iostream"
#include <GL/glu.h>
#include <GL/glut.h>
#include <datastruct.h>
class Window:public QOpenGLWindow, protected QOpenGLFunctions
{
    Q_OBJECT
public:
    ~Window();
    void initializeGL();
    void resizeGL(int width,int height);
    void paintGL();
    void keyPressEvent(QKeyEvent *e);
    void mousePressEvent(QMouseEvent *f);
    void mouseMoveEvent(QMouseEvent *g);
    void mouseReleaseEvent(QMouseEvent *h);
    void wheelEvent(QWheelEvent *l);
    //void readfile();
    void DrawFlight();
    int numpoint;
    GLfloat spanx,spany,spanz;
    GLfloat Scalef;
    mesh* flight;
    QPointF tempP,tempR;
    bool pressflag;
    bool fullscreen;
    GLfloat  triangle_rotate ;
    GLfloat  quads_rotate ;
    GLfloat  circle_rotate ;
    QVector3D eye,center,up;
};

#endif // WINDOW_H
