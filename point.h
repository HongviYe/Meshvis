#ifndef POINT_H
#define POINT_H
#include <QtOpenGL>

class Point
{
public:
    Point(GLfloat x1,GLfloat y1,GLfloat z1);
    Point();
    void DarwLine(Point P);
//private:
    GLfloat x,y,z;
};

#endif // POINT_H
