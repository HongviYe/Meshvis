#include<fstream>
#include<map>
#include<iostream>
#include<algorithm>
#include<cmath>
#include <assert.h>
#include<string>
#include"geom_func.h"
#include"datastruct.h"
using namespace std;
double Plane::Aera() {
	double s[3],a=0;
	for (int i = 0; i < 3; i++) {
		s[i] = Distance(*p[(i + 1) % 3], *p[(i + 2) % 3]);
		a += s[i];
	}
	a /= 2;
	return sqrt(a*(a-s[0])*(a-s[1])*(a-s[2]));


}
/*when the tone vtk point coordinate do not match the other's, we need it to store one parentcell of every point*/
mesh::~mesh() {
	delete[] Point_array;
	delete[] Tetra_array;
	delete[] Plane_array;
	if(beita_array)
	delete[] beita_array;
	if(temp_times_array)
	delete[] temp_times_array;
	PlanePoint_map.clear();
}
void Boundaryline::OpenGLDrawLine(){
    glBegin(GL_LINES);
    glVertex3f(p[0]->coordinate[0],p[0]->coordinate[1],p[0]->coordinate[2]);
    glVertex3f(p[1]->coordinate[0],p[1]->coordinate[1],p[1]->coordinate[2]);
    glEnd();
}
int mesh::readVTK(const char *fname1, const char *fname2) {
	ifstream finsurface(fname2);
	ifstream finspace(fname1);
    assert(finspace);
    assert(finsurface);
	string Trash;
	Point P;
    Boundaryline bl;
	int temp, index[4];
	Plane S;
	PlanePoint SP;
	Tetrahedron T;
    while (Trash != "POINTS"){
		finspace >> Trash;
        cout<<Trash;
    }
	finspace >> numpoint;
	Point_array = new Boundarypoint[numpoint];
	finspace >> Trash;

	/*get point_array*/

	for (int i = 0; i < numpoint; i++) {
		for (int j = 0; j < 3; j++)
			finspace >> Point_array[i].coordinate[j];
	}
	finspace >> Trash;
	finspace >> numtetra;
	Tetra_array = new  Tetrahedron[numtetra];
	assert(Tetra_array);
	Plane_array = new Plane[4*numtetra];
	assert(Plane_array);
	finspace >> Trash;

	/*get tetrahedron and the plane on it*/
	for (int i = 0; i < numtetra; i++) {
		finspace >> temp;
		for (int j = 0; j < 4; j++) {
			finspace >> index[j];
            T.p[j] = &Point_array[index[j]];
		}
		Tetra_array[i] = T;
		for (int j = 0; j < 4; j++) {
			SP.p[0] = index[(j + 1) % 4];
			SP.p[1] = index[(j + 2) % 4];
			SP.p[2] = index[(j + 3) % 4];
			sort(SP.p, SP.p + 3);	
			
			/*if the plane had not read yet*/
			if (PlanePoint_map.find(SP) == PlanePoint_map.end()) {
				S.p[0] = &Point_array[SP.p[0]];
				S.p[1] = &Point_array[SP.p[1]];
				S.p[2] = &Point_array[SP.p[2]];
				S.Left = &Tetra_array[i];
				S.Right = NULL;
				Plane_array[numsurface++] = S;
/*
                for(int t=0;t<3;t++){
                    int k;
                    for(k=0;k<S.p[t]->nei_B_Point.size();k++){
                        if(S.p[(t)%3]->nei_B_Point[k]==S.p[(t+1)%3])
                            break;
                    }
                    if(k==S.p[(t)%3]->nei_B_Point.size()){
                        S.p[(t+1)%3]->nei_B_Point.push_back(S.p[(t)%3]);
                        S.p[(t)%3]->nei_B_Point.push_back(S.p[(t+1)%3]);
                        bl.set(S.p[(t+1)%3],S.p[(t)%3]);
                        BL.push_back(bl);
                    }

                }
*/
				PlanePoint_map.insert(pair<PlanePoint, Plane*>(SP, &Plane_array[numsurface-1]));
			}
			else { /*if the plane had read already*/
				PlanePoint_map[SP]->Right = &Tetra_array[i];
			}
			Tetra_array[i].s[j] = PlanePoint_map[SP];
			Tetra_array[i].color = 0;
		}

/*
        for(int j=0;j<3;j++){
            int k;
            for(k=0;k<Point_array[index[(j+1)%3]].nei_B_Point.size();k++){
                if(Point_array[index[(j+1)%4]].nei_B_Point[k]==&Point_array[index[(j)%4]])
                    break;
            }
            if(k==Point_array[index[(j+1)%3]].nei_B_Point.size()){
                Point_array[index[(j+1)%4]].nei_B_Point.push_back(&Point_array[index[(j)%4]]);
                Point_array[index[(j)%4]].nei_B_Point.push_back(&Point_array[index[(j+1)%4]]);
                bl.set(&Point_array[index[(j+1)%4]],&Point_array[index[(j)%4]]);
                BL.push_back(bl);
            }

        }
*/    }

	finspace.close();
	/* end reading of spacemesh */
	while (Trash != "POINTS") {
		finsurface >> Trash;
	}
	finsurface >> numPointBoundary;
	
	finsurface >> Trash;
	for (int i = 0; i < numPointBoundary; i++) {
		for (int j = 0; j < 3; j++) {
			finsurface >> Trash;
			P.coordinate[j] = Point_array[i].coordinate[j];
		}
	}
	finsurface >> Trash;
	finsurface >> numPlaneBoundary;
	Planeboundary_array = new Plane[numPlaneBoundary];
	finsurface >> Trash;
	for (int i = 0; i <  numPlaneBoundary; i++) {
		finsurface >> Trash;
		int te;
		for (int k = 0; k < 3; k++) {
			finsurface>>te;
			S.p[k] = &Point_array[te];
		}

		Planeboundary_array[i] = S;
        for(int t=0;t<3;t++){
            int k;
            for(k=0;k<S.p[t]->nei_B_Point.size();k++){
                if(S.p[(t)%3]->nei_B_Point[k]==S.p[(t+1)%3])
                    break;
            }
            if(k==S.p[(t)%3]->nei_B_Point.size()){
                S.p[(t+1)%3]->nei_B_Point.push_back(S.p[(t)%3]);
                S.p[(t)%3]->nei_B_Point.push_back(S.p[(t+1)%3]);
                bl.set(S.p[(t+1)%3],S.p[(t)%3]);
                BL.push_back(bl);
            }

        }

	}
	while (Trash != "float") 
		finsurface >> Trash;
	for (int i = 0; i < numPointBoundary; i++) {
		for (int j = 0; j < 3; j++) {
			finsurface >> BP[i].direction[j];
		}
	}
	/*
	for (int i = 0; i < numtetra; i++) {
		for (int j = 0; j < 4; j++) {
			if (PointBoundaryPoint_map.find(*Tetra_array[i].p[j]) != PointBoundaryPoint_map.end()) {
				PointBoundaryPoint_map[*Tetra_array[i].p[j]]->oneparentcell = &Tetra_array[i];
			}
		}
	}
	*/


	/*delete the map and close the file stream obj*/
	finsurface.close();
	return 0;
}
inline const double Orient3d(const Point P,const Plane S) {
	double temp=0;

	/* if the point (X,Y,Z,A) has axis (a,b,c), A is the top of tetradron
	*  we do X=X-A,Y=Y-A,Z=Z-A,then the volumn of tetradron is: 
	*         | Xa  Ya  Za |
	*         |            |
	* V= 1/36 | Xb  Yb  Zb |
	*         |            |
	*         | Xc  Yc  Zc |
	*/
	Point PP[3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)
			PP[i].coordinate[j]=S.p[i]->coordinate[j] - P.coordinate[j];
	}
	temp = (PP[0].coordinate[0]* PP[1].coordinate[1]* PP[2].coordinate[2]) 
		+ (PP[0].coordinate[1] * PP[1].coordinate[2] * PP[2].coordinate[0]) 
		+ (PP[0].coordinate[2] * PP[1].coordinate[0] * PP[2].coordinate[1])
		- (PP[0].coordinate[2] * PP[1].coordinate[1] * PP[2].coordinate[0])
		- (PP[1].coordinate[0] * PP[0].coordinate[1] * PP[2].coordinate[2])
		- (PP[0].coordinate[0] * PP[2].coordinate[1] * PP[1].coordinate[2]);

	/*recover the data*/
	return 1.0/36*temp;
}
/*inline function to reduce the cost of calling function*/
inline double Distance(const Point P1, const Point P2) {
	double temp=0;
	for (int i = 0; i < 3; i++)
		temp += (P1.coordinate[i] - P2.coordinate[i])*(P1.coordinate[i] - P2.coordinate[i]);
	return sqrt(temp);
}
//#ifdef _DEBUG
double ObjectFunction(Boundarypoint root_P,Point remote_P,Point P,Plane S) {
	double ans = Distance(P, remote_P);
	double temp0 = 0, temp1 = 0, temp2, temp3;
	for (int i = 0; i < 3; i++) {
		temp0 += (P.coordinate[i]-root_P.coordinate[i]) * root_P.direction[i];
		temp1 += (P.coordinate[i] - root_P.coordinate[i])*(P.coordinate[i] - root_P.coordinate[i]);
	}
	temp0 *= temp0;
	temp3 = abs(Orient3d(P,S));
	temp2 = exp(sqrt(abs(temp1-temp0)+1)/5);
	if (temp2 >1e8)
		temp2 = 1e8;
	ans +=temp2;
	ans += 1.0 / pow(temp3, 1);
	return ans;
}
/*the below function just for debug*/
void mesh::FindBadPoint() {
	double thick = 0.1*(1-pow(1.1,30))/(1-1.1);
	ofstream fout1("output.vtk", ios::app);
	fout1 << "POINT_DATA " << numPointBoundary << endl;
	fout1 << "SCALARS aveH int" << endl;
	fout1 << "LOOKUP_TABLE H" << endl;
	for (int i = 0; i < numPointBoundary; i++) {
		if (0.9*BP[i].MMD > 2 * thick )
			fout1 << 0 << endl;
		else
			fout1 << 1 << endl;
	}
	fout1.close();
}
void mesh::aveH() {
	temp_times_array=new int[numPointBoundary];
	int index[3];
	for (int i = 0; i < numPointBoundary; i++) {
		BP[i].H = 0;
		temp_times_array[i] = 0;
	}
	ifstream finsurface("surfacemesh.vtk");
	string s;
	while (s != "CELLS")
		finsurface >> s;
	finsurface >> s;
	finsurface >> s;
	for (int i = 0; i < numPlaneBoundary; i++) {
		finsurface >> s;
		for (int j = 0; j < 3; j++)
			finsurface >> index[j];
		for (int j = 0; j < 3; j++) {
			BP[index[(j + 1) % 3]].H += Distance(BP[index[j % 3]], BP[index[(j+1) % 3]]);
			BP[index[(j ) % 3]].H += Distance(BP[index[j % 3]], BP[index[(j + 1) % 3]]);
			temp_times_array[index[(j+1) % 3]]++;
			temp_times_array[index[j % 3]]++;
		}
	}
	for (int i = 0; i < numPointBoundary; i++) {
		BP[i].H = BP[i].H / temp_times_array[i];
		if (BP[i].H > 1.25*BP[i].MMD*0.1) {
			BP[i].H = 1.25*BP[i].MMD*0.1;
		}
	}
}
//#endif
vector<Tetrahedron*>  findshell(Point* lastPointA,Point* lastPointB, Tetrahedron* lastTetradron,bool isboundary) {
	vector<Tetrahedron*> ans;
	isboundary = false;
	int fir = 0;
	Tetrahedron* cur1=lastTetradron,*cur2=lastTetradron;
	while (lastTetradron != cur1||fir==0) {
		fir = 1;
		for (int i = 0; i < 4; i++) {
			int c = 0;
			for (int j = 0; j < 3; j++) {
				if (lastTetradron->s[i]->p[j] == lastPointA || lastTetradron->s[i]->p[j] == lastPointB) {
					c++;
				}
			}
			if (c == 2) {
				int flag = 1;
				if (lastTetradron->s[i]->Left != lastTetradron&&lastTetradron->s[i]->Left != cur2) {
					ans.push_back(lastTetradron->s[i]->Left);
					cur2 = lastTetradron;
					lastTetradron = lastTetradron->s[i]->Left;
					flag = 0;
				}
				else if (lastTetradron->s[i]->Right != NULL&&lastTetradron->s[i]->Right != cur2) {
					ans.push_back(lastTetradron->s[i]->Right);
					cur2 = lastTetradron;
					lastTetradron = lastTetradron->s[i]->Right;
					flag = 0;
				}
				else if (lastTetradron->s[i]->Right == NULL)
					isboundary = true;
				if(flag)
				break;
			}
		}
	}
	return ans;
}
int mesh::MMDslove1() {
	Point remotePoint;
	double volumn = 0;
	Plane *currentplane,*bestplane,*lastsurface;
	Tetrahedron *currentTetradron,*bestTetradron,*lastTetradron;
	double mindistance;
	int i;
#pragma omp parallel for num_threads(4) private(currentplane,bestplane,lastsurface,currentTetradron,bestTetradron,lastTetradron,mindistance,remotePoint,volumn)
	for (i = 0; i <numpoint; i++) {
		int c = 1;
		/**
		for (int k = 0; k < 4; k++) {
			onep[i]->s[k]->Left->color = 2;
			if (onep[i]->s[k]->Right != NULL)
				onep[i]->s[k]->Right->color = 2;
			else
				cout << 222;
		}//*/
		if(BP[i].direction[1]== BP[i].direction[0]&& BP[i].direction[2] == BP[i].direction[0]&& BP[i].direction[0]==0)
			continue;
		for (int j = 0; j < 3; j++)
				remotePoint.coordinate[j] = BP[i].coordinate[j] + MAXSIZE*BP[i].direction[j];
		//lastTetradron = BP[i].oneparentcell;
		lastTetradron = BP[i].onefather;
		if (lastTetradron == NULL) {
			cerr << "there is one cell do not have a parentcell"<<endl;
			continue;
		}
		double ave = 1.0/100*abs(Orient3d(*lastTetradron->p[0],*lastTetradron->s[0]));
		currentTetradron = NULL;
		lastsurface = NULL;
		currentplane = NULL;
		while (1) {
			bestplane = NULL; bestTetradron = NULL;/*to store the best choice of the jumping*/
			int k;
			lastTetradron->color= c++;
			mindistance = inf;
			volumn = 0;
			for (int j = 0; j < 4; j++) {
				if (lastTetradron->s[j] != lastsurface) {
					if (Orient3d(*lastTetradron->p[j], *lastTetradron->s[j])*Orient3d(remotePoint, *lastTetradron->s[j]) <= 0) {
						currentplane = lastTetradron->s[j];
						/*to get the neighber cell*/
						if (currentplane->Left != lastTetradron) {
							currentTetradron = currentplane->Left;
						}
						else if (currentplane->Right != NULL) {
							currentTetradron = currentplane->Right;
						}
						else {
							if (volumn < abs(Orient3d(Point(BP[i]), *currentplane)))
								volumn = abs(Orient3d(Point(BP[i]), *currentplane));
							continue;
						}
						/*hybrid object funtion*/
						if (abs(Orient3d(Point(BP[i]), *currentplane)) >= volumn) {

						for (k = 0; k < 3; k++) {
							if (currentTetradron->s[k] == currentplane)
								break;
						}
						if (ObjectFunction(BP[i], remotePoint, *currentTetradron->p[k],*currentTetradron->s[k]) < mindistance) {
							mindistance = ObjectFunction(BP[i], remotePoint, *currentTetradron->p[k], *currentTetradron->s[k]);
							bestplane = currentplane;
							bestTetradron = currentTetradron;
						}
						}
					}
				}
			}
			/*if we can not find plane fit the need*/
			if (bestplane == NULL) {
				double linep[2][3], facep[3][3], intPnt[3];
				int intType, intCod;
				bool bEpsilon=false;
				Point PP;
				mindistance = inf;
				bool flag = true;
				for (int k = 0; k < 4; k++) {
					if (lastTetradron->s[k]->Right == NULL) {
						for (int t = 0; t < 3; t++) {
							for (int y = 0; y < 3; y++)
								facep[t][y] = lastTetradron->s[k]->p[t]->coordinate[y];
						}
						for (int k = 0; k < 3; k++) {
							linep[0][k] = BP[i].coordinate[k];
							linep[1][k] = BP[i].coordinate[k] + MAXSIZE * BP[i].direction[k];
						}
						GEOM_FUNC::lin_tri_intersect3d(linep, facep, &intType, &intCod, intPnt, bEpsilon);
						for (int y = 0; y < 3; y++)
							PP.coordinate[y] = intPnt[y];
						if (intType == GEOM_FUNC::LTI_INTERSECT_FAC) {
							BP[i].MMD = Distance(BP[i], PP);
							flag = false;
						}
						for (int j = 0; j < 3; j++) {
							if (Distance(BP[i], *lastTetradron->s[k]->p[j]) < mindistance)
								mindistance = Distance(BP[i], *lastTetradron->s[k]->p[j]);
						}
					}
				}
				if (flag)
					BP[i].MMD = mindistance;
				break;
			}
			/*if we find at lest one corresponding plane*/
			lastsurface = bestplane;
			lastTetradron = bestTetradron;
			for (k = 0; k < 3; k++) {
				if (bestTetradron->s[k] == bestplane) {
					break;
				}
			}
			//if(PointBoundaryPoint_map.find(*bestTetradron->p[k])!= PointBoundaryPoint_map.end())
			//PointBoundaryPoint_map[*bestTetradron->p[k]]->MMD = c;
		}
	}
	return 0;
}
int mesh::MMDslove2() {
	bool bEpsilon=false,flag=false;
	bool isread[4] = {false};
	int lastsituation ;
	Point remotePoint,P,*lastPointA,*lastPointB;
	Plane *lastsurface;
	Tetrahedron *lastTetradron;
	double linep[2][3], facep[3][3], intPnt[3];
	int intType, intCod;
	for (int i = 0; i < numtetra; i++) {
		for (int t = 0; t < 4; t++) {
			for (int j = 0; j < 3; j++) {
				linep[0][j] = Tetra_array[i].p[t]->coordinate[j];
				linep[1][j] = Tetra_array[i].p[t]->coordinate[j] + MAXSIZE*Tetra_array[i].p[t]->direction[j];
				for (int k = 0; k < 3; k++) {
					facep[j][k] = Tetra_array[i].s[t]->p[j]->coordinate[k];
				}
			}
			GEOM_FUNC::lin_tri_intersect3d(linep, facep, &intType, &intCod, intPnt, bEpsilon);
			if (intType ==GEOM_FUNC::LTI_INTERSECT_FAC) {
				Tetra_array[i].p[t]->onefather = &Tetra_array[i];
			}
		}
	}
	
	
	for (int i = 0; i < numPointBoundary; i++) {
		flag = false;
		lastsurface = NULL;
		lastsituation = -1;
		lastTetradron = BP[i].onefather;
		if (lastTetradron == NULL)
			continue;
		for (int j = 0; j < 3; j++)
			remotePoint.coordinate[j] = BP[i].coordinate[j] + MAXSIZE*BP[i].direction[j];
		for (int j = 0; j < 3; j++) {
			linep[0][j] = BP[i].coordinate[j];
			linep[1][j] = BP[i].coordinate[j] + MAXSIZE*BP[i].direction[j];
		}
		while (true) {
			bool isboundary;
			int nodecount = 0;
			int edgecount = 0;
			for (int k = 0; k < 4; k++) {
				for (int m = 0; m < 3; m++)
					for (int n = 0; n < 3; n++)
						facep[m][n] = lastTetradron->s[k]->p[m]->coordinate[n];
				if (lastTetradron->s[k] != lastsurface) {
					GEOM_FUNC::lin_tri_intersect3d(linep, facep, &intType, &intCod, intPnt, bEpsilon);
					if (intType == GEOM_FUNC::LTI_INTERSECT_FAC) {
						if (lastTetradron->s[k]->Right == NULL) {
							for (int t = 0; t < 3; t++)
								P.coordinate[t] = intPnt[t];
							BP[i].MMD = Distance(P, BP[i]);
							flag = true;
						}
						else {
							lastsurface = lastTetradron->s[k];
							if (lastTetradron->s[k]->Left != lastTetradron)
								lastTetradron = lastTetradron->s[k]->Left;
							else
								lastTetradron = lastTetradron->s[k]->Right;
							lastsituation = GEOM_FUNC::LTI_INTERSECT_FAC;
							lastPointA = NULL;
							lastPointB = NULL;
						}
					}
					else if (intType == GEOM_FUNC::LTI_INTERSECT_EDG) {
						edgecount++; 
						Point *C, *D,*E;
						isread[k] = true;
						if (edgecount == 2 ) {
							if (lastsituation == GEOM_FUNC::LTI_INTERSECT_EDG) {
								for (int m = 0; m < 4; m++) {
									if (isread[m] == true) {
										int c = 0;
										for (int n = 0; n < 3; n++) {
											if (lastTetradron->s[m]->p[n] == lastPointA || lastTetradron->s[m]->p[n] == lastPointB) {
												c++;
												E = lastTetradron->s[m]->p[n];
											}
										}
										if (c == 2)
											C = lastTetradron->p[m];
										if (c == 1)
											D = E;
										
									}
								}
								lastPointA = C;
								lastPointB = D;
							}
							else{
								int c = 0;
								for (int m = 0; m < 4; m++) {
									if (isread[m] == false) {
										c++;
										if (c == 1)
											lastPointA = lastTetradron->p[m];
										else
											lastPointB = lastTetradron->p[m];
									}

								}
							}
							vector<Tetrahedron*> shell = findshell(lastPointA, lastPointB , lastTetradron,isboundary);
							if (isboundary) {
								for (int t = 0; t < 3; t++)
									P.coordinate[t] = intPnt[t];
								BP[i].MMD = Distance(P, BP[i]);
								flag = true;
							}
							else {
								for (unsigned int t = 0; t < shell.size(); t++) {
									for (int y = 0; y < 3; y++) {
										for (int m = 0; m < 3; m++)
											for (int n = 0; n < 3; n++) 
												facep[m][n] = shell[t]->s[y]->p[m]->coordinate[n];
										GEOM_FUNC::lin_tri_intersect3d(linep, facep, &intType, &intCod, intPnt, bEpsilon);
										if (intType == GEOM_FUNC::LTI_INTERSECT_INS || intType == GEOM_FUNC::LTI_INTERSECT_FAC){
											lastTetradron = shell[t];
										}
										if (lastTetradron == shell[t] && intType == GEOM_FUNC::LTI_INTERSECT_EDG)
											lastsurface = shell[t]->s[y];

									}
								}
							}
							lastsituation = GEOM_FUNC::LTI_INTERSECT_EDG;
						}
					}
					else if (intType==GEOM_FUNC::LTI_INTERSECT_NOD) {
						nodecount++;
						if (nodecount == 2&&lastsurface!=NULL) {
							for (int t = 0; t < 3; t++)
								P.coordinate[t] = intPnt[t];
							BP[i].MMD = Distance(P, BP[i]);
							flag = true;
						}
					}
				}
			}
			if (flag)
				break;
		}
	}
	for (int i = 0; i < numPointBoundary; i++)
		if (BP[i].MMD == 0)
			BP[i].MMD = 20000;
	return 0;
}
int mesh::writeVTK(const char* fname2) {
	/*ouput mmd into file*/
	ifstream fin(fname2);
	ofstream fout("output.vtk");
	char buffer[256];
	while (!fin.eof()) {
		fin.read(buffer,256);
		int n = fin.gcount();
		fout.write(buffer,n);
	}
	fout.close();
	ofstream fout1("output.vtk",ios::app);
	fout1 << "SCALARS mmd float" << endl;
	fout1 << "LOOKUP_TABLE MMD"<<endl;
	for (int i = 0; i < numPointBoundary; i++)
		fout1 << (abs(BP[i].MMD)) << endl;
	if (beita_array) {
		fout1 << "CELL_DATA " << numPlaneBoundary << endl << "SCALARS beita float 1" << endl;
		fout1 << "LOOKUP_TABLE default2" << endl;
		for (int i = 0; i < numPlaneBoundary; i++)
			fout1 << beita_array[i] << endl;
	}
	/*
	fout1 << "CELL_DATA " << numtetra << endl << "SCALARS celltype float 1" << endl;
	fout1 << "LOOKUP_TABLE default1" << endl;
	for (int i = 0; i < numtetra; i++)
		fout1 << Tetra_array[i].color << endl;
		//*/
	return 0;
}
void mesh::FindBeta() {
	beita_array = new double[numPlaneBoundary];
	for (int i = 0; i < numPlaneBoundary; i++)
		beita_array[i] = 0;
	double H[3],temp[3],ans;
	double K[3][3];
	for (int i = 0; i < numPlaneBoundary; i++) {
		ans = 0;
		for (int j =0 ; j < 3; j++) {
			H[j] = Planeboundary_array[i].p[j]->H;
			temp[j] = 0;
			for (int k = 0; k < 3; k++) {
				K[j][k] = 0;
			}
		}
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int t = 0; t < 3; t++) {
					K[j][k] += (Planeboundary_array[i].p[(j + 1) % 3]->coordinate[t] - Planeboundary_array[i].p[(j + 2) % 3]->coordinate[t]) * (Planeboundary_array[i].p[(k + 1) % 3]->coordinate[t] - Planeboundary_array[i].p[(k + 2) % 3]->coordinate[t]);
				}
			}
		}
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				K[j][k] /= 4 * pow(Planeboundary_array[i].Aera(),2);
			}
		}
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				temp[j] += H[k] * K[j][k];
			}
		}
		for (int j = 0; j < 3; j++)
			ans += temp[j] * H[j];
		beita_array[i] = exp(sqrt(ans+1e-13));
	}
}
#ifdef _DEBUG
int mesh::accurate() {
	for (int i = 0; i < numPointBoundary; i++)
		BP[i].MMD = 20000;
	Point temp;
	bool bEpsilon=false;
	bool flag;
	int i;
#pragma omp parallel for num_threads(4) private(temp,bEpsilon,flag)
	for (int i = 0; i <numPointBoundary; i++) {
		if(i<numPointBoundary/4-1&&i%100==0)
		cout << i*100.0 / numPointBoundary << endl;
		double linep[2][3], facep[3][3], intPnt[3];
		int intType, intCod;
		double mindistance = inf;
		for (int k = 0; k < 3; k++) {
			linep[0][k] = BP[i].coordinate[k];
			linep[1][k] = BP[i].coordinate[k] + 40000*BP[i].direction[k];
		}
		
		for (int j = 0; j < numPlaneBoundary; j++) {
			flag = false;
			for (int k = 0; k < 3; k++) {
				if(*Planeboundary_array[j].p[k] == Point_array[i]){
					flag = true;
					break;
				}
			}
			if (flag==true)
				continue;
			for (int t = 0; t < 3; t++) {
				for (int k = 0; k < 3; k++)
					facep[t][k] = Planeboundary_array[j].p[t]->coordinate[k];
			}
			GEOM_FUNC::lin_tri_intersect3d(linep, facep, &intType, &intCod, intPnt, bEpsilon);
			if (intType == GEOM_FUNC::LTI_INTERSECT_FAC) {
				for (int k = 0; k < 3; k++)
					temp.coordinate[k] = intPnt[k];
				if (BP[i].MMD > Distance(BP[i], temp))
				BP[i].MMD = Distance(BP[i],temp);
			}
		}
	}
	return 0;
}
#endif
