//
// Created by paola on 06/11/19.
//

#include "MyLibrary.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include </usr/local/MATLAB/R2019b/extern/include/mat.h>
#include </usr/local/MATLAB/R2019b/extern/include/matrix.h>

#define pi M_PI
void readFile(double* xa, double* xb, double* ya, double* yb, double* za, double* zb, double* h, double* tini, double* tfin, double* dt)
{
    FILE *fp;
    char stringa80[80]; //Stringa per lettura file
    fp = fopen("Files/input.txt", "r");
    fgets(stringa80, 80, fp);
    sscanf(stringa80, "%lf %lf", xa, xb);
    fgets(stringa80, 80, fp);
    sscanf(stringa80, "%lf %lf", ya, yb);
    fgets(stringa80, 80, fp);
    sscanf(stringa80, "%lf %lf", za, zb);
    fgets(stringa80, 80, fp);
    sscanf(stringa80, "%lf ",h);
    fgets(stringa80, 80, fp);
    sscanf(stringa80, "%lf %lf %lf",tini,tfin,dt);
    fclose(fp);

}

void input_broadcast(double* xa, double* xb, double* ya, double* yb, double* za, double* zb, double* h, double* tini, double* tfin, double* dt, int root)
{
    MPI_Bcast(xa,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD);
    MPI_Bcast(xb,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD);
    MPI_Bcast(ya, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD);
    MPI_Bcast(yb, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD);
    MPI_Bcast(za, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD);
    MPI_Bcast(zb, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD);
    MPI_Bcast(h,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD);
    MPI_Bcast(tini,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD);
    MPI_Bcast(tfin,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD);
    MPI_Bcast(dt,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD);

}

void position(int myid, int cuts, int* xpos,int* ypos,int* zpos)
{
    int cuts2 = cuts*cuts;
    *xpos = myid%cuts;
    *ypos = (int) (myid % cuts2)/cuts;
    *zpos = (int) myid/cuts2;
}

void findNeighbours(int myid,int cuts,int xpos,int ypos,int zpos ,int* xneigh1,int* xneigh2,int* yneigh1,int* yneigh2,int* zneigh1,int* zneigh2)
{
    *xneigh1 = myid -1;
    *xneigh2 = myid +1;
    *yneigh1 = myid - cuts;
    *yneigh2 = myid + cuts;
    *zneigh1 = myid - (cuts*cuts);
    *zneigh2 = myid + (cuts*cuts);

    if(xpos == 0) *xneigh1 = -1;
    if(xpos == cuts-1) *xneigh2 = -1;
    if(ypos == 0) *yneigh1 = -1;
    if(ypos == cuts-1) *yneigh2 = -1;
    if(zpos == 0) *zneigh1 = -1;
    if(zpos == cuts-1) *zneigh2 = -1;
}

void cubecoord(int cuts, double h, int xpos, int ypos, int zpos, double xa, double xb, double ya, double yb,
               double za, double zb, int* nx, int* ny, int* nz)
{
    *nx = (int) (fabs(xa-xb)/h)+1;
    *ny = (int) (fabs(ya-yb)/h)+1;
    *nz = (int) (fabs(za-zb)/h)+1;

    cubeloc.nxloc = *nx/cuts;
    cubeloc.nyloc = *ny/cuts;
    cubeloc.nzloc = *nz/cuts;

    cubeloc.x0 = xa + (xpos *  cubeloc.nxloc) * h;
    cubeloc.x1 = xa + (xpos *  cubeloc.nxloc +  cubeloc.nxloc - 1) * h;
    cubeloc.y0 = ya + (ypos *  cubeloc.nyloc) * h;
    cubeloc.y1 = ya + (ypos *  cubeloc.nyloc +  cubeloc.nyloc - 1) * h;
    cubeloc.z0 = za + (zpos *  cubeloc.nzloc) * h;
    cubeloc.z1 = za + (zpos *  cubeloc.nzloc +  cubeloc.nzloc - 1) * h;

    if(xpos==(cuts-1))
    {
        cubeloc.nxloc = *nx-cubeloc.nxloc*(cuts-1);
        cubeloc.x1 = xb;
    }
    if(ypos==(cuts-1))
    {
        cubeloc.nyloc = *ny-cubeloc.nyloc*(cuts-1);
        cubeloc.y1 = yb;
    }
    if(zpos==(cuts-1))
    {
        cubeloc.nzloc = *nz-cubeloc.nzloc*(cuts-1);
        cubeloc.z1 = zb;
    }
}

void allocations(double**** u0, double**** u1,Cubeloc cubeloc){
    double ***u1bis = NULL, ***u0bis = NULL;
    int i = 0, j = 0;

    u0bis = (double***) malloc(sizeof(double**)*cubeloc.nxloc);
    for(i = 0; i < cubeloc.nxloc; i++){
        u0bis[i] = (double**) malloc(sizeof(double*)*cubeloc.nyloc);
        for(j = 0; j < cubeloc.nyloc; j++){
            u0bis[i][j] = (double*)calloc(cubeloc.nzloc,sizeof(double));
        }
    }
    *u0 = u0bis;

    u1bis = (double***) malloc(sizeof(double**)*cubeloc.nxloc);
    for(i = 0; i < cubeloc.nxloc; i++){
        u1bis[i] = (double**) malloc(sizeof(double*)*cubeloc.nyloc);
        for(j = 0; j < cubeloc.nyloc; j++){
            u1bis[i][j] = (double*)calloc(cubeloc.nzloc,sizeof(double));
        }
    }
    *u1 = u1bis;

}

double exSolution(double x, double y, double z, double t)
{
    return sin(2*pi*(x+y+z-t));
}

double f(double x,double y, double z,double t)
{
    return -2*pi*cos(2*pi*(x+y+z-t))+12*pi*pi*sin(2*pi*(x+y+z-t));
}

void exportToParaview(double u1[cubeloc.nxloc][cubeloc.nyloc][cubeloc.nzloc])
{
    int i=0,j=0,k=0;
    FILE *file;
    //fp=freopen(NULL,"w+",fp);

    file = fopen("Files/Output.csv","w");
    fprintf(file,"x coord, y coord, z coord, scalar\n");
    for (i = 0; i < cubeloc.nxloc; i++) {
        for (j = 0; j < cubeloc.nyloc; j++) {
            for (k = 0; k < cubeloc.nzloc; k++) {
                fprintf(file, "%d,%d,%d,%lf\n", i, j, k, u1[i][j][k]);
            }
        }
    }
    fclose(file);
}

void MatlabFilesClean(){
    int myid = 0;
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    MATFile *pmat;

    char id[12];
    sprintf(id, "%d", myid);

    char filename[80];
    strcpy(filename, "Files/MatlabScreens/id_");
    strcat(filename, id);
    strcat(filename,".mat");

    pmat = matOpen(filename, "w");
    matClose(pmat);
}

void MatlabScreen(double u[cubeloc.nxloc][cubeloc.nyloc][cubeloc.nzloc], int its){
    int myid = 0;
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    MATFile *pmat;
    mxArray *pa1;
    const mwSize dims[3] = {(mwSize)cubeloc.nzloc,(mwSize)cubeloc.nyloc,(mwSize)cubeloc.nxloc};

    char id[12];
    sprintf(id, "%d", myid);
    char iter[12];
    sprintf(iter, "%d", its);

    char filename[80];
    strcpy(filename, "Files/MatlabScreens/id_");
    strcat(filename, id);
    strcat(filename,".mat");

    char variablename[80];
    strcpy(variablename, "u2_");
    strcat(variablename, id);
    /*
    strcat(variablename, "_iter_");
    strcat(variablename, iter);
*/
    pmat = matOpen(filename, "u");

    pa1 = mxCreateNumericArray((mwSize)3,dims,mxDOUBLE_CLASS,mxREAL);
    memcpy((void *)(mxGetPr(pa1)), (void *)u, sizeof(double)*(cubeloc.nxloc*cubeloc.nyloc*cubeloc.nzloc));
    matPutVariable(pmat,variablename,pa1);

    matClose(pmat);
}

void sendRecvAllocation(double*** uSend,double*** uRecv,int x,int y)
{
    int i=0;
    double **uSend2 = (double**)malloc(sizeof(double*)*x);
    for(i=0;i<x;i++)
    {
        uSend2[i]=(double*)calloc(y,sizeof(double));
    }
    *uSend=uSend2;

    double **uRecv2 = (double**)malloc(sizeof(double*)*x);
    for(i=0;i<x;i++)
    {
        uRecv2[i]=(double*)calloc(y,sizeof(double));
    }
    *uRecv=uRecv2;
}

void stampaArray(double u1[cubeloc.nxloc][cubeloc.nyloc][cubeloc.nzloc])
{
    int i=0,j=0,k=0;
    FILE *file;
    //fp=freopen(NULL,"w+",fp);

    file = fopen("Files/Stampe.csv","w");
    fprintf(file,"x coord, y coord, z coord, scalar\n");
    for (i = 0; i < cubeloc.nxloc; i++) {
        for (j = 0; j < cubeloc.nyloc; j++) {
            for (k = 0; k < cubeloc.nzloc; k++) {
                fprintf(file, "%d,%d,%d,%lf\n", i, j, k, u1[i][j][k]);
            }
        }
    }
    fclose(file);
}