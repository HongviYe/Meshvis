#include "Window.h"
#include <QtGui>
#include <QtCore>
#include <QtOpenGL>
#include <QTimer>
#include <QTransform>
#include <fstream>
#include <string>
#include <cmath>
#define GL_PI 3.1415926
#define GL_RADIUX  0.2f
Window::~Window()
{
        triangle_rotate = 0.0;
        quads_rotate = 0.0;
        circle_rotate = 0.0;
}
void Window::DrawFlight(){
    GLfloat cof=1e5;
    for(int i=0;i<flight->BL.size();i++){
        glBegin(GL_LINES);
        glVertex3f(flight->BL[i].p[0]->coordinate[0]/cof,flight->BL[i].p[0]->coordinate[1]/cof,flight->BL[i].p[0]->coordinate[2]/cof);
        glVertex3f(flight->BL[i].p[1]->coordinate[0]/cof,flight->BL[i].p[1]->coordinate[1]/cof,flight->BL[i].p[1]->coordinate[2]/cof);
        glEnd();
    }
}
void Window::initializeGL(){
    initializeOpenGLFunctions();
    setGeometry(300, 150, 500, 500);//设置窗口初始位置和大小
        glShadeModel(GL_SMOOTH);//设置阴影平滑模式
        glClearColor(0.0, 0.0, 0.0, 0.0);//改变窗口的背景颜色，不过我这里貌似设置后并没有什么效果
        glClearDepth(1.0);//设置深度缓存    for(int i=0;i<flight->BL.size();i++){

        glEnable(GL_DEPTH_TEST);//允许深度测试
        glDepthFunc(GL_LEQUAL);//设置深度测试类型
        glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);//进行透视校正
        spanx=spany=spanz=10;
        Scalef=1;pressflag=false;
        up.setX(1.0);
        up.setY(0.0);
        up.setZ(0.0);
        eye.setX(0.3);
        eye.setY(0.3);
        eye.setZ(0.3);
        center.setX(0.3);
        center.setY(0.0);
        center.setZ(0.0);
}
void Window::resizeGL(int width, int height){
    resize(QSize(width,height));
    glMatrixMode(GL_PROJECTION);//选择投影矩阵
       glLoadIdentity();//重置选择好的投影矩阵
       //gluPerspective(90.0, (GLfloat)width/(GLfloat)height, 0.1, 100.0);//建立透视投影矩阵
       gluLookAt(eye.x(),eye.y(),eye.z(),center.x(),center.y(),center.z(),up.x(),up.y(),up.z());//center_poi.x(),center_poi.y(),0);
       glMatrixMode(GL_MODELVIEW);
       glLoadIdentity();
}
void Window::paintGL(){

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //glClear()函数在这里就是对initializeGL()函数                                                   //中设置的颜色和缓存深度等起作用
    glLoadIdentity();//重置当前的模型观察矩阵,该句执行完后，将焦点移动到了屏幕的中心
    //gluLookAt(10,10,10,0,0,0,1,1,1);
    gluLookAt(eye.x(),eye.y(),eye.z(),center.x(),center.y(),center.z(),up.x(),up.y(),up.z());
    glRotatef(spanx,0.0,1.0,0.0);
    glRotatef(spany,1.0,0.0,0.0);//*/
    glScalef(Scalef,Scalef,Scalef);
    /*
    QVector3D eye(22222.0,22222.0,22222.0),center(2,2,0),up(0.3,0.8,-0.5);
    QMatrix4x4 transform;
    transform.setToIdentity();gluLookAt(eye.x(),eye.y(),eye.z(),center.x(),center.y(),center.z(),up.x(),up.y(),up.z());
    transform.lookAt(eye,center,up);//*/
    DrawFlight();


}
 void Window::keyPressEvent(QKeyEvent *e){
     switch(e->key())
         {
             //F1键为全屏和普通屏显示切换键
             case Qt::Key_F1:
                 fullscreen = !fullscreen;
                 if(fullscreen)
                     showFullScreen();
                 else
                 {
                     showNormal();
                 }
                 //updateGL();
                 break;
             //Ese为退出程序键
             case Qt::Key_Escape:
                 close();


         }
 }
void Window::mousePressEvent(QMouseEvent *e){
    if(e->button()==Qt::LeftButton){
        pressflag=true;
        tempP=e->screenPos();
        update();
    }
    if(e->button()==Qt::MidButton){

        update();
    }

}
void Window::wheelEvent(QWheelEvent *e){
        Scalef+=e->delta()/1e3;
        update();
}
void Window::mouseReleaseEvent(QMouseEvent *h){
    if(h->button()==Qt::LeftButton){
        pressflag=false;
       tempR=h->screenPos();
         update();
    }
    if(h->button()==Qt::RightButton){
        Scalef=1;
        update();
    }
}
void Window::mouseMoveEvent(QMouseEvent *g){
    GLfloat sca;
    if(g->isAccepted()){
        if(pressflag){
           eye.setX(eye.x()+sin(g->pos().rx()/800));
            spany=-g->pos().ry();

        }
    }

    update();
}


/*
void Window::readfile(){

    //glLineWidth(0.3);//设置线段宽度
    double maxsize=1;
    std::string trash("i");
    std::ifstream fin("/home/yhf/test/sample");
    if(!fin){
        std::cout<<"can not open the file;";
        return;
    }

    while(trash!="POINTS")
        fin>>trash;
    fin>>numpoint;
    P=new Point[numpoint];
    for(int i=0;i<numpoint;i++){
        fin>>P[i].x>>P[i].y>>P[i].z;
        if(abs(P[i].x)>maxsize)
            maxsize=abs(P[i].x);
    }
    maxsize=12000;
    for(int i=0;i<numpoint;i++){
        P[i].x/=maxsize;
        P[i].y/=maxsize;
        P[i].z/=maxsize;
    }
    fin>>trash;
    fin>>numline;
    fin>>trash;
    linearray=new int[3*numline];
    //fin>>trash;
    for(int i=0;i<numline;i++){
        fin>>trash;
        for(int j=0;j<3;j++)
               fin>>linearray[3*i+j];
    }
}
*/
