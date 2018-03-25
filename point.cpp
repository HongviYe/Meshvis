#include "point.h"

Point::Point(GLfloat x1,GLfloat y1,GLfloat z1)
{
        x=x1;y=y1;z=z1;
}
Point::Point()
{
        x=y=z=0;

}
void Point::DarwLine(Point P){
      glBegin(GL_LINES);
      glColor3f(1.0,0.0,0.0);
      glVertex3f(x,y,z); //定点坐标范围
      glVertex3f(P.x,P.y,P.z);
      glEnd();
}
