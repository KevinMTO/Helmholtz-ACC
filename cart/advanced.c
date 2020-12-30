/* 
**************************************************************************
*                                                                        *
* Subroutine to solve the Helmholtz equation:                            *
* (d2/dx2)u + (d2/dy2)u - alpha u = f                                    *
*                                                                        *
* Solves poisson equation on rectangular grid assuming:                  *
* (1) Uniform discretization in each direction, and                      *
* (2) Dirichlect boundary conditions.                                    *
* Jacobi iterative method is used in this routine.                       *
*                                                                        *
* Input:  n,m         Number of grid points in the X/Y directions        *
*         dx,dy       Grid spacing in the X/Y directions                 *
*         alpha       Helmholtz eqn. coefficient                         *
*         omega       Relaxation factor                                  *
*         f(n,m)      Right hand side function                           *
*         u(n,m)      Dependent variable (solution)                      *
*         tolerance   Tolerance for iterative solver                     *
*         maxit       Maximum number of iterations                       *
*                                                                        *
* Output: u(n,m) - Solution                                              *
*                                                                        *
**************************************************************************
*/

#include <mpi.h>
#if defined(WIN32) || defined (_CLUSTER_OPENMP)
/* fix a visualstudio 2005 bug */
#include <omp.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "jacobi_enhance.h"

#define U(j,i) afU[((j) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((i)- data->iColFirst)]
#define F(j,i) afF[((j) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((i)- data->iColFirst)]
#define UOLD(j,i) uold[ ((j) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((i)- data->iColFirst)]

extern void ExchangeJacobiMpiData(struct JacobiData *data, double *uold, double *afU) ;

void Jacobi (struct JacobiData *data)
{
    /*use local pointers for performance reasons*/
    double *afU, *afF, *temp;
    int i, j;
    double fLRes;
    
    double ax, ay, b, residual, tmpResd;
    
    double *uold = (double*) calloc(
        ( (data->iColLast - data->iColFirst +1) * (data->iRowLast - data->iRowFirst + 1) ), sizeof(double));
    afU = data->afU;
    afF = data->afF;


    /*
    #pragma omp parallel for private(j, i) 
    for (j = data->iRowFirst ; j <= data->iRowLast ; j++)
    {
        for (i = data->iColFirst ; i <= data->iColLast ; i++)
        {
            UOLD(j,i)=0.0;
        }
    }

    if(data->iMyRank==0){
    printf("\n000 %d %d %d %d %d %d\n",data->iRowFirst, data->iColFirst,data->iRowLast, data->iColLast,(data->iRowLast - data->iRowFirst),(data->iColLast- data->iColFirst));
    UOLD((data->iRowFirst), data->iColFirst)=1.0;
    printf("%d ind0ice1\n",( (((data->iRowFirst)) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((data->iColFirst)- data->iColFirst)));
    UOLD((data->iRowFirst), data->iColLast)=1.0;
    printf("%d ind0ice2\n",( (((data->iRowFirst)) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((data->iColLast)- data->iColFirst)));
    UOLD((data->iRowLast), data->iColFirst)=1.0;
    printf("%d ind0ice3\n",( (((data->iRowLast)) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((data->iColFirst)- data->iColFirst)));
    UOLD((data->iRowLast), data->iColLast)=1.0;
    printf("%d ind0ice4\n",( (((data->iRowLast)) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((data->iColLast)- data->iColFirst)));
    }
    if(data->iMyRank==2){
    printf("\n222 %d %d %d %d %d %d\n",data->iRowFirst, data->iColFirst,data->iRowLast, data->iColLast,(data->iRowLast - data->iRowFirst),(data->iColLast- data->iColFirst));
    UOLD((data->iRowFirst), data->iColFirst)=1.0;
    printf("%d ind2ice1\n",( (((data->iRowFirst)) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((data->iColFirst)- data->iColFirst)));
    UOLD((data->iRowFirst), data->iColLast)=1.0;
    printf("%d ind2ice2\n",( (((data->iRowFirst)) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((data->iColLast)- data->iColFirst)));
    UOLD((data->iRowLast), data->iColFirst)=1.0;
    printf("%d ind2ice3\n",( (((data->iRowLast)) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((data->iColFirst)- data->iColFirst)));
    UOLD((data->iRowLast), data->iColLast)=1.0;
    printf("%d ind2ice4\n",( (((data->iRowLast)) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + ((data->iColLast)- data->iColFirst)));
    }
    //check if iniztialed
    if(data->iMyRank==0){
        printf("UOLD OF 0\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n ", uold[0], uold[1], uold[2], uold[3], uold[4],uold[5],uold[6],uold[7],uold[8],uold[9],uold[10],uold[11],uold[12],uold[13],uold[14],uold[15],uold[16],uold[17],uold[18],uold[19],uold[20],uold[21],uold[22],uold[23],uold[24],uold[25],uold[26],uold[27],uold[28],uold[29],uold[30],uold[31],uold[32],uold[33],uold[34],uold[35]);
    } 
    if(data->iMyRank==2){
        printf("UOLD OF 2\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n ", uold[0], uold[1], uold[2], uold[3], uold[4],uold[5],uold[6],uold[7],uold[8],uold[9],uold[10],uold[11],uold[12],uold[13],uold[14],uold[15],uold[16],uold[17],uold[18],uold[19],uold[20],uold[21],uold[22],uold[23],uold[24],uold[25],uold[26],uold[27],uold[28],uold[29],uold[30],uold[31],uold[32],uold[33],uold[34],uold[35]);
    } 
    if(data->iMyRank==3){
        printf("UOLD OF 3\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n ", uold[0], uold[1], uold[2], uold[3], uold[4],uold[5],uold[6],uold[7],uold[8],uold[9],uold[10],uold[11],uold[12],uold[13],uold[14],uold[15],uold[16],uold[17],uold[18],uold[19],uold[20],uold[21],uold[22],uold[23],uold[24],uold[25],uold[26],uold[27],uold[28],uold[29],uold[30],uold[31],uold[32],uold[33],uold[34],uold[35]);
    } 
        //check if iniztialed
    if(data->iMyRank==0){
        printf("afU OF 0\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n ", afU[0], afU[1], afU[2], afU[3], afU[4],afU[5],afU[6],afU[7],afU[8],afU[9],afU[10],afU[11],afU[12],afU[13],afU[14],afU[15],afU[16],afU[17],afU[18],afU[19],afU[20],afU[21],afU[22],afU[23],afU[24],afU[25],afU[26],afU[27],afU[28],afU[29],afU[30],afU[31],afU[32],afU[33],afU[34],afU[35]);
    } 
    if(data->iMyRank==2){
        printf("afU OF 2\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n ", afU[0], afU[1], afU[2], afU[3], afU[4],afU[5],afU[6],afU[7],afU[8],afU[9],afU[10],afU[11],afU[12],afU[13],afU[14],afU[15],afU[16],afU[17],afU[18],afU[19],afU[20],afU[21],afU[22],afU[23],afU[24],afU[25],afU[26],afU[27],afU[28],afU[29],afU[30],afU[31],afU[32],afU[33],afU[34],afU[35]);
    } 
    if(data->iMyRank==3){
        printf("afU OF 3\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n ", afU[0], afU[1], afU[2], afU[3], afU[4],afU[5],afU[6],afU[7],afU[8],afU[9],afU[10],afU[11],afU[12],afU[13],afU[14],afU[15],afU[16],afU[17],afU[18],afU[19],afU[20],afU[21],afU[22],afU[23],afU[24],afU[25],afU[26],afU[27],afU[28],afU[29],afU[30],afU[31],afU[32],afU[33],afU[34],afU[35]);
    } 
    //to remove
    int how_many=0;

    */

    if (uold)
    {
        ax = 1.0 / (data->fDx * data->fDx);      /* X-direction coef */
        ay = 1.0 / (data->fDy * data->fDy);      /* Y_direction coef */
        b = -2.0 * (ax + ay) - data->fAlpha;     /* Central coeff */
        residual = 10.0 * data->fTolerance;

        while (data->iIterCount < data->iIterMax && residual > data->fTolerance) 
        {
            residual = 0.0;
            
            temp = uold;
            uold = afU;
            afU = temp;
            

            ExchangeJacobiMpiData(data, uold, afU); 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
        if(how_many<6){
            if(data->iMyRank==0){
                printf("CHECK MATRIX 0-1\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %d\n ", uold[0], uold[1], uold[2], uold[3], uold[4],uold[5],uold[6],uold[7],uold[8],uold[9],uold[10],uold[11],uold[12],uold[13],uold[14],uold[15],uold[16],uold[17],uold[18],uold[19],uold[20],uold[21],uold[22],uold[23],uold[24],uold[25],uold[26],uold[27],uold[28],uold[29],uold[30],uold[31],uold[32],uold[33],uold[34],uold[35],how_many);
            } 
            if(data->iMyRank==2){
                printf("CHECK MATRIX 2-1\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %d\n ", uold[0], uold[1], uold[2], uold[3], uold[4],uold[5],uold[6],uold[7],uold[8],uold[9],uold[10],uold[11],uold[12],uold[13],uold[14],uold[15],uold[16],uold[17],uold[18],uold[19],uold[20],uold[21],uold[22],uold[23],uold[24],uold[25],uold[26],uold[27],uold[28],uold[29],uold[30],uold[31],uold[32],uold[33],uold[34],uold[35],how_many);
            } 
            if(data->iMyRank==3){
                printf("CHECK MATRIX 3-1\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %d\n", uold[0], uold[1], uold[2], uold[3], uold[4],uold[5],uold[6],uold[7],uold[8],uold[9],uold[10],uold[11],uold[12],uold[13],uold[14],uold[15],uold[16],uold[17],uold[18],uold[19],uold[20],uold[21],uold[22],uold[23],uold[24],uold[25],uold[26],uold[27],uold[28],uold[29],uold[30],uold[31],uold[32],uold[33],uold[34],uold[35],how_many);
            } 
        }
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	        #pragma omp parallel for private(j, i, fLRes) reduction(+:residual)
            for (j = data->iRowFirst + 1; j <= data->iRowLast - 1; j++)
            {
                for (i = data->iColFirst + 1; i <= data->iColLast - 1; i++)
                {
                    fLRes = ( ax * (UOLD(j, i-1) + UOLD(j, i+1))
                            + ay * (UOLD(j-1, i) + UOLD(j+1, i))
                            +  b * UOLD(j, i) - F(j, i)) / b;

                    //UPDATE SOLUTION
                    U(j,i) = UOLD(j,i) - data->fRelax * fLRes;
                    
                    /* accumulate residual error */
                    residual += fLRes * fLRes;
                }
            }
            /* end omp for */
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
        if(how_many<6){
            if(data->iMyRank==0){
                printf("CHECK MATRIX 0-2\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %d\n ", uold[0], uold[1], uold[2], uold[3], uold[4],uold[5],uold[6],uold[7],uold[8],uold[9],uold[10],uold[11],uold[12],uold[13],uold[14],uold[15],uold[16],uold[17],uold[18],uold[19],uold[20],uold[21],uold[22],uold[23],uold[24],uold[25],uold[26],uold[27],uold[28],uold[29],uold[30],uold[31],uold[32],uold[33],uold[34],uold[35],how_many);
            } 
            if(data->iMyRank==2){
                printf("CHECK MATRIX 2-2\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %d\n ", uold[0], uold[1], uold[2], uold[3], uold[4],uold[5],uold[6],uold[7],uold[8],uold[9],uold[10],uold[11],uold[12],uold[13],uold[14],uold[15],uold[16],uold[17],uold[18],uold[19],uold[20],uold[21],uold[22],uold[23],uold[24],uold[25],uold[26],uold[27],uold[28],uold[29],uold[30],uold[31],uold[32],uold[33],uold[34],uold[35],how_many);
            } 
            if(data->iMyRank==3){
                printf("CHECK MATRIX 3-2\n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %lf  %lf  %lf  %lf  %lf  %lf  \n %d\n ", uold[0], uold[1], uold[2], uold[3], uold[4],uold[5],uold[6],uold[7],uold[8],uold[9],uold[10],uold[11],uold[12],uold[13],uold[14],uold[15],uold[16],uold[17],uold[18],uold[19],uold[20],uold[21],uold[22],uold[23],uold[24],uold[25],uold[26],uold[27],uold[28],uold[29],uold[30],uold[31],uold[32],uold[33],uold[34],uold[35],how_many);
            } 
        
            how_many++;
        }
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            tmpResd = residual;
            MPI_Allreduce(&tmpResd, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            /* error check */
            (data->iIterCount)++;
            residual = sqrt(residual) / (data->iCols * data->iRows);

            

            
        }
        //while

        data->fResidual = residual;
        free(uold);
    }
    else 
    {
        fprintf(stderr,"Error: can't allocate memory\n");
        Finish(data);
        exit(1);
    }
}

void ExchangeJacobiMpiData(struct JacobiData *data, double *uold, double * afU){
    //printf("MPI %d CAME BY my proc NULL is %d\n\n\n", data->iMyRank, MPI_PROC_NULL);
    

    //assign communicator for topology
    MPI_Comm topocomm = data->Topo_comm;
    //printf("MPI: %d TOPO COMM %lu\n",data->iMyRank,topocomm);

    MPI_Request request[8];
    MPI_Status  status[8];
    
    double *afF;
    int iReqCnt = 0;
    int i, j;

    const int iTagMoveDown = 8;
    const int iTagMoveUp = 9;
    const int iTagMoveWest = 10;
    const int iTagMoveEast = 11;

    
    afF = data->afF;

    int north, west, south, east;
    
    MPI_Cart_shift(topocomm, 0, 1, &west, &east);
    MPI_Cart_shift(topocomm, 1, 1, &north, &south);
    
    //printf("MPI: %d has north %d south %d east %d west %d\n",data->iMyRank,north,south,east,west);

    int width = (data->iColLast - data->iColFirst + 1); 
    int height = (data->iRowLast - data->iRowFirst + 1);


    
////////TYPES
    int n_blocks = height-2, stride = width; // -2 because I want to avoid race conditions on corners
    MPI_Datatype column_type;

    //Create a derived type
    MPI_Type_vector(n_blocks, 1, stride, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);

    MPI_Datatype row_type;

    //Create a derived type
    MPI_Type_contiguous((width-2), MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);


    if(west != MPI_PROC_NULL){


        // Send a column
        //printf("MPI %d Swest\n",data->iMyRank);

        MPI_Isend(&UOLD((data->iRowFirst+1), (data->iColFirst + 1)), 1, column_type, west, iTagMoveWest, topocomm, &request[iReqCnt]);
  
        iReqCnt++;
        //printf("MPI %d west\n",data->iMyRank);
        MPI_Irecv(&UOLD((data->iRowFirst+1), data->iColFirst), 1, column_type, west, iTagMoveEast, topocomm,&request[iReqCnt]);
        iReqCnt++;


    }

    if(east != MPI_PROC_NULL){



        // Send a column
	    //printf("MPI %d Seast\n",data->iMyRank);

        MPI_Isend(&UOLD((data->iRowFirst+1), (data->iColLast-1)), 1, column_type, east, iTagMoveEast, topocomm, &request[iReqCnt]);
        iReqCnt++;
        //printf("MPI %d east\n",data->iMyRank);

        MPI_Irecv(&UOLD((data->iRowFirst+1), data->iColLast), 1, column_type, east, iTagMoveWest, topocomm,&request[iReqCnt]);
        iReqCnt++;


    }

    if(north != MPI_PROC_NULL){
	        
        //printf("MPI %d Snorth\n",data->iMyRank);

        MPI_Isend(&UOLD( (data->iRowFirst + 1), (data->iColFirst + 1)), 1, row_type, north, iTagMoveUp, topocomm, &request[iReqCnt]);
        iReqCnt++;
        //printf("MPI %d north\n",data->iMyRank);

	    MPI_Irecv(&UOLD( (data->iRowFirst), (data->iColFirst + 1)), 1, row_type, north, iTagMoveDown, topocomm, &request[iReqCnt]);
        iReqCnt++;	

    }
    if(south != MPI_PROC_NULL){
        
        //printf("MPI %d Ssouth\n",data->iMyRank);

        MPI_Isend(&UOLD( (data->iRowLast - 1), (data->iColFirst+1)), 1, row_type, south, iTagMoveDown, topocomm, &request[iReqCnt]);
        iReqCnt++;
        //printf("MPI %d south\n",data->iMyRank);

	    MPI_Irecv(&UOLD( (data->iRowLast), (data->iColFirst+1)), 1, row_type, south, iTagMoveUp, topocomm, &request[iReqCnt]);
        iReqCnt++;
    }

    /*
    #pragma omp parallel for private(i)
    for (j = data->iRowFirst + 1; j <= data->iRowLast - 1; j++)
    {
        for (i = data->iColFirst + 1; i <= data->iColLast - 1; i++)
        {
            UOLD(j, i) = U(j, i);
        }
    }
    */
    

    MPI_Waitall(iReqCnt, request, status);
    //Free the type after it is not needed
    MPI_Type_free(&column_type);
    MPI_Type_free(&row_type);
    //printf("MPI %d EXITS FROM EXCHANGE\n",data->iMyRank);
}

