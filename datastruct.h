#pragma once
#include<vector>
#include<cmath>
#include<iostream>
#include<map>
#include <QtOpenGL>
#ifndef _mmd_datastruct_H
#define _mmd_datastruct_H
using namespace std;
class Tetrahedron;
/*to store the point obj*/
class Point {
public:
	double coordinate[3];
	bool operator < (const Point SP)const {
		for (int i = 0; i<3; i++)
			if (SP.coordinate[i] != coordinate[i]) {
				return SP.coordinate[i] > coordinate[i];
			}
		return false;
	}
	bool operator > (const Point SP)const {
		return SP < *this;
	}
	bool operator==(const Point SP)const {
		return !(SP > *this || SP < *this);
	}
};
class Boundarypoint:public Point {
public:
	double direction[3];
	//Tetrahedron* oneparentcell;
	double MMD;
	double H;
    vector<Boundarypoint*> nei_B_Point;
	Tetrahedron* onefather;
};

/*to hash the point_index to plane*/
class PlanePoint {
public:
	int p[3];
	bool operator < (const PlanePoint SP)const {
		for(int i=0;i<3;i++)
		if (SP.p[i] != p[i]) {
			return SP.p[i] > p[i];
		}
		return false;
	}
	bool operator > (const PlanePoint SP)const{
		return SP < *this;
	}
	bool operator==(const PlanePoint SP)const {
		return !(SP > *this || SP < *this);
	}

};
/*to store the plane obj*/
struct Plane {
	Boundarypoint* p[3];
	Tetrahedron* Left;
	Tetrahedron* Right;
	double Aera();
};
class Boundaryline{
public:
    Boundaryline(){

    }
    Boundaryline(Boundarypoint* p1,Boundarypoint* p0){
       p[0]=p0;
       p[1]=p1;
    }
    void set(Boundarypoint* p1,Boundarypoint* p0){
       p[0]=p0;
       p[1]=p1;
    }
    Boundarypoint* p[2];
    void OpenGLDrawLine();
};
/*to store Tetrahedron obj*/
class Tetrahedron {
public:
	Plane* s[4];
	Boundarypoint* p[4];
	double color;
};
class mesh {
public:
	const double MAXSIZE ;
	const double inf ;
	Boundarypoint *BP;
	int *temp_times_array;
	/*when the tone vtk point coordinate do not match the other's, we need it to store one parentcell of every point*/
	Tetrahedron *Tetra_array;
	Plane *Plane_array;
    Boundarypoint* &Point_array ;
	Plane *Planeboundary_array;
	map <PlanePoint, Plane*> PlanePoint_map;
    int numpoint, numsurface , numtetra, &numPointBoundary, numPlaneBoundary;
    double *beita_array;
    vector<Boundaryline> BL;
	/**user input**/
	double ThickOfFirstLevel;
	double retio;
    double Numlevel;

	/**************/
public:
    mesh(const char *fname1, const char *fname2):MAXSIZE(100000.0),inf(1e20),beita_array(NULL),Point_array(BP),numsurface(0), numPointBoundary (numpoint) {
		readVTK(fname1,fname2);
	}
	~mesh();
	int readVTK(const char *fname1, const char *fname2);
	int MMDslove1();
	int MMDslove2();

	int writeVTK(const char* fname2);
	int accurate();
	void aveH();
	void FindBadPoint();
	void FindBeta();
};
extern inline double Distance(const Point P1, const Point P2);
#endif // !_HEAD
