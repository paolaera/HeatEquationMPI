#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "cmake-build-debug/MyLibrary.h"

#define CACCA printf("Dall'id %d ho ricevuto cacca\n",myid);

Cubeloc cubeloc;

int main()
{
    //Dichiarazione delle variabili
    double xa = 0, xb = 0, ya = 0, yb = 0, za = 0, zb = 0;
    double h = 0,tini = 0,tfin = 0,dt = 0,nts=0;
    int myid = -1,np = 0,root = 0;
    int nx = 0, ny =0, nz = 0,cuts = 0;
    int xneigh1=0,xneigh2=0,yneigh1=0,yneigh2=0,zneigh1=0,zneigh2=0,xpos=0,ypos=0,zpos=0,nxloc=0,nyloc=0,nzloc=0;
    double **xSend1 = NULL, **xSend2 = NULL, **ySend1 = NULL, **ySend2 = NULL, **zSend1 = NULL, **zSend2 = NULL;
    double **xRecv1 = NULL, **xRecv2 = NULL, **yRecv1 = NULL, **yRecv2 = NULL, **zRecv1 = NULL, **zRecv2 = NULL;
    int i=0, j=0 , k=0,its=0;
    double effe = 0, t = 0, h2=0;

    //Inizializzazione di MPI
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&np);

    MatlabFilesClean();
    //Lettura del file contenente la dimensione della grigliau
    if(myid==root)
    {
        readFile(&xa, &xb, &ya, &yb, &za, &zb, &h, &tini, &tfin, &dt);
    }

    input_broadcast(&xa, &xb, &ya, &yb, &za, &zb, &h, &tini, &tfin, &dt, root);

    //printf("%lf \t %lf in myid=%d\n",x0,x1,myid);


    //Numero di punti della discretizzazione
    nx = (int)(fabs(xb-xa)/h)+1;
    ny = (int)(fabs(yb - ya) / h) + 1;        //+1 per far tornare il conteggio dei punti
    nz = (int)(fabs(zb - za) / h) + 1;
    cuts =(int) pow(np,1.00/3);
    /*if(myid==root)
    {
        printf("%d \t %d",np,cuts);
    }*/

    //Ogni processo sa la sua posizione e sa con chi confina
    position(myid,cuts,&xpos,&ypos,&zpos);
    findNeighbours(myid,cuts, xpos, ypos, zpos, &xneigh1, &xneigh2, &yneigh1, &yneigh2, &zneigh1, &zneigh2);
    //Riempimento dei minicubi
    cubecoord(cuts,h,xpos,ypos,zpos,xa,xb,ya,yb,za,zb,&nx,&ny,&nz);

    //Creazione degli array per le iterazioni
    double u0[cubeloc.nxloc][cubeloc.nyloc][cubeloc.nzloc];
    double u1[cubeloc.nxloc][cubeloc.nyloc][cubeloc.nzloc];

    for(i=0;i<cubeloc.nxloc;i++)
    {
        for(j=0;j<cubeloc.nyloc;j++)
        {
            for(k=0;k<cubeloc.nzloc;k++)
            {

                u0[i][j][k] =exSolution(cubeloc.x0+i*h,cubeloc.y0+j*h,cubeloc.z0+k*h,tini);
                u1[i][j][k] =exSolution(cubeloc.x0+i*h,cubeloc.y0+j*h,cubeloc.z0+k*h,tini+dt);

            }
        }
    }



    //allocazione matrici send e recv
    sendRecvAllocation(&xSend1,&xRecv1,cubeloc.nyloc,cubeloc.nzloc);
    sendRecvAllocation(&xSend2,&xRecv2,cubeloc.nyloc,cubeloc.nzloc);

    sendRecvAllocation(&ySend1,&yRecv1,cubeloc.nzloc,cubeloc.nyloc);
    sendRecvAllocation(&ySend2,&yRecv2,cubeloc.nzloc,cubeloc.nyloc);

    sendRecvAllocation(&zSend1,&zRecv1,cubeloc.nxloc,cubeloc.nyloc);
    sendRecvAllocation(&zSend2,&zRecv2,cubeloc.nxloc,cubeloc.nyloc);

    //RiempimentoMi son

    nts =(int)((tfin-tini)/dt);
    h2=h*h;
    for(its=0;its<nts;its++) {
        t = tini + its * dt;
        if (xneigh1 != -1) {
            //Riempio xSend1
            for (i = 0; i < cubeloc.nyloc; i++) {
                for (j = 0; j < cubeloc.nzloc; j++) {
                    xSend1[i][j] = u1[0][i][j];
                }
                //Mando xSend1
                MPI_Send(xSend1[i], cubeloc.nzloc, MPI_DOUBLE, xneigh1, i, MPI_COMM_WORLD);
            }
            //Ricevo xRecv1
            for (i = 0; i < cubeloc.nyloc; i++) {
                MPI_Recv(xRecv1[i], cubeloc.nzloc, MPI_DOUBLE, xneigh1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        if (xneigh2 != -1) {
            //Riempio xSend2
            for (i = 0; i < cubeloc.nyloc; i++) {
                for (j = 0; j < cubeloc.nzloc; j++) {
                    xSend2[i][j] = u1[cubeloc.nxloc - 1][i][j];
                }
                //Mando xSend2
                MPI_Send(xSend2[i], cubeloc.nzloc, MPI_DOUBLE, xneigh2, i, MPI_COMM_WORLD);
            }
            //Ricevo xRecv2
            for (i = 0; i < cubeloc.nyloc; i++) {
                MPI_Recv(xRecv2[i], cubeloc.nzloc, MPI_DOUBLE, xneigh2, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        if (yneigh1 != -1) {
            //Riempio ySend1
            for (i = 0; i < cubeloc.nxloc; i++) {
                for (j = 0; j < cubeloc.nzloc; j++) {
                    ySend1[i][j] = u1[i][0][j];
                }
                //Mando ySend1
                MPI_Send(ySend1[i], cubeloc.nzloc, MPI_DOUBLE, yneigh1, i, MPI_COMM_WORLD);
            }
            //Ricevo yRecv1
            for (i = 0; i < cubeloc.nxloc; i++) {
                MPI_Recv(yRecv1[i], cubeloc.nzloc, MPI_DOUBLE, yneigh1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        if (yneigh2 != -1) {
            //Riempio ySend2
            for (i = 0; i < cubeloc.nxloc; i++) {
                for (j = 0; j < cubeloc.nzloc; j++) {
                    ySend2[i][j] = u1[i][cubeloc.nyloc - 1][j];
                }
                //Mando ySend2
                MPI_Send(ySend2[i], cubeloc.nzloc, MPI_DOUBLE, yneigh2, i, MPI_COMM_WORLD);
            }
            //Ricevo xRecv2
            for (i = 0; i < cubeloc.nxloc; i++) {
                MPI_Recv(yRecv2[i], cubeloc.nzloc, MPI_DOUBLE, yneigh2, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        if (zneigh1 != -1) {
            //Riempio zSend1
            for (i = 0; i < cubeloc.nxloc; i++) {
                for (j = 0; j < cubeloc.nyloc; j++) {
                    zSend1[i][j] = u1[i][j][0];
                }
                //Mando zSend1
                MPI_Send(zSend1[i], cubeloc.nyloc, MPI_DOUBLE, zneigh1, i, MPI_COMM_WORLD);
            }
            //Ricevo zRecv1
            for (i = 0; i < cubeloc.nxloc; i++) {
                MPI_Recv(zRecv1[i], cubeloc.nyloc, MPI_DOUBLE, zneigh1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        if (zneigh2 != -1) {
            //Riempio zSend2
            for (i = 0; i < cubeloc.nxloc; i++) {
                for (j = 0; j < cubeloc.nyloc; j++) {
                    zSend2[i][j] = u1[i][j][cubeloc.nzloc - 1];
                }
                //Mando zSend2
                MPI_Send(zSend2[i], cubeloc.nyloc, MPI_DOUBLE, zneigh2, i, MPI_COMM_WORLD);
            }
            //Ricevo zRecv2
            for (i = 0; i < cubeloc.nxloc; i++) {
                MPI_Recv(zRecv2[i], cubeloc.nyloc, MPI_DOUBLE, zneigh2, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        for (i = 1; i < cubeloc.nxloc - 1; i++) {
            for (j = 1; j < cubeloc.nyloc - 1; j++) {
                for (k = 1; k < cubeloc.nzloc - 1; k++) {

                    effe = f(cubeloc.x0 + i * h, cubeloc.y0 + j * h, cubeloc.z0 + k * h, t);
                    u1[i][j][k] = u0[i][j][k] + dt * (effe + 1 / h2 * (-6 * u0[i][j][k]));

                    if (i != 0) u1[i][j][k] += (dt / h2) * u0[i - 1][j][k];
                    else u1[i][j][k] += (dt / h2) * xRecv1[j][k];

                    if (i != cubeloc.nxloc - 1) u1[i][j][k] += (dt / h2) * u0[i + 1][j][k];
                    else u1[i][j][k] += (dt / h2) * xRecv2[j][k];

                    if (j != 0) u1[i][j][k] += (dt / h2) * u0[i][j - 1][k];
                    else u1[i][j][k] += (dt / h2) * yRecv1[i][k];

                    if (j != cubeloc.nyloc - 1) u1[i][j][k] += (dt / h2) * u0[i][j + 1][k];
                    else u1[i][j][k] += (dt / h2) * yRecv2[i][k];

                    if (k != 0) u1[i][j][k] += (dt / h2) * u0[i][j][k - 1];
                    else u1[i][j][k] += (dt / h2) * zRecv1[i][j];

                    if (k != cubeloc.nzloc - 1) u1[i][j][k] += (dt / h2) * u0[i][j][k + 1];
                    else u1[i][j][k] += (dt / h2) * zRecv2[i][j];

                    if (i == 0) {
                        if (xneigh1 == -1)
                            u1[0][j][k] = exSolution(cubeloc.x0, cubeloc.y0 + j * h, cubeloc.z0 + k * h, t);
                        if (xneigh2 == -1)
                            u1[cubeloc.nxloc - 1][j][k] = exSolution(cubeloc.x1, cubeloc.y0 + j * h, cubeloc.z0 + k * h,
                                                                     t);
                    }

                    if (j == 0) {
                        if (yneigh1 == -1)
                            u1[i][0][k] = exSolution(cubeloc.x0 + i * h, cubeloc.y0, cubeloc.z0 + k * h, t);
                        if (yneigh2 == -1)
                            u1[i][cubeloc.nyloc - 1][k] = exSolution(cubeloc.x0 + i * h, cubeloc.y1, cubeloc.z0 + k * h,
                                                                     t);
                    }

                    if (k == 0) {
                        if (zneigh1 == -1)
                            u1[i][j][0] = exSolution(cubeloc.x0 + i * h, cubeloc.y0 + j * h, cubeloc.z0, t);
                        if (zneigh2 == -1)
                            u1[i][j][cubeloc.nzloc - 1] = exSolution(cubeloc.x0 + i * h, cubeloc.y0 + j * h, cubeloc.z1, t);
                    }
                }
            }
            printf("Dall'id=%d ho ottenuto u2[0]"
        }
        for(i=0;i<cubeloc.nxloc;i++){
            for(j=0;j<cubeloc.nyloc;j++){
                for(k=0;k<cubeloc.nzloc;k++){
                    u0[i][j][k] = u1[i][j][k];
                }
            }
        }
    }
    exportToParaview(u1);

/*
    if(myid==7) {
        for (i = 0; i < cubeloc.nxloc; i++) {
            for (j = 0; j < cubeloc.nyloc; j++) {
                for (k = 0; k < cubeloc.nzloc; k++) {
                    printf("u0[%d][%d][%d]=%lf\n", i, j, k, u0[i][j][k]);
                }
            }
        }
    }
*/

    MPI_Finalize();
    return 0;

}

