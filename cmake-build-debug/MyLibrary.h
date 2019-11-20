//
// Created by paola on 06/11/19.
//

#ifndef HEATEQUATION_MYLIBRARY_H
#define HEATEQUATION_MYLIBRARY_H

#endif //HEATEQUATION_MYLIBRARY_H

typedef struct {
    double x0;
    double x1;
    double y0;
    double y1;
    double z0;
    double z1;
    int nxloc;
    int nyloc;
    int nzloc;
}Cubeloc;

extern Cubeloc cubeloc;
void readFile(double* xa, double* xb, double* ya, double* yb, double* za, double* zb, double* h, double* tini, double* tfin, double* dt);
void input_broadcast(double* xa, double* xb, double* ya, double* yb, double* za, double* zb, double* h, double* tini, double* tfin, double* dt, int root);
void position(int myid,int cuts,int* xpos,int* ypos,int* zpos);
void findNeighbours(int myid,int cuts,int xpos,int ypos,int zpos ,int* xneigh1,int* xneigh2,int* yneigh1,int* yneigh2,int* zneigh1,int* zneigh2);
void cubecoord(int cuts, double h, int xpos, int ypos, int zpos, double xa, double xb, double ya, double yb,
               double za, double zb, int* nx, int* ny, int* nz);
void allocations(double**** u0, double**** u1,Cubeloc cubeloc);
double exSolution(double x, double y, double z, double t);
double f(double x,double y, double z,double t);
void exportToParaview(double u0[cubeloc.nxloc][cubeloc.nyloc][cubeloc.nzloc]);
void MatlabScreen(double u[cubeloc.nxloc][cubeloc.nyloc][cubeloc.nzloc], int its);
void MatlabFilesClean();
void sendRecvAllocation(double*** uSend,double*** uRecv,int x,int y);


