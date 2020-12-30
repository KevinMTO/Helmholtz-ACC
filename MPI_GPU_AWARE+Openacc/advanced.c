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

#include <openacc.h>
#include <mpi-ext.h>
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

#define LENGTH ( (data->iColLast - data->iColFirst +1) * (data->iRowLast - data->iRowFirst + 1) )

extern void ExchangeJacobiMpiData(struct JacobiData *data, double *  uold, double *  afU) ;

void set_local_rank(struct JacobiData * restrict data, int * local);







void Jacobi (struct JacobiData * restrict data)
{
    /*use local pointers for performance reasons*/
    double * restrict afU;
    double * restrict afF;
    double * restrict temp;

    int i, j;
    double fLRes;

    int localRank = -1;


    
    double ax, ay, b, residual, tmpResd, residualup, residualdown, residualwest, residualeast, residualinner;
    
    double * restrict uold =  calloc(
        ( (data->iColLast - data->iColFirst +1) * (data->iRowLast - data->iRowFirst + 1) ), sizeof(double));

    afU = data->afU;
    afF = data->afF;


    if (uold)
    {

        ax = 1.0 / (data->fDx * data->fDx);      /* X-direction coef */
        ay = 1.0 / (data->fDy * data->fDy);      /* Y_direction coef */
        b = -2.0 * (ax + ay) - data->fAlpha;     /* Central coeff */
        residual = 10.0 * data->fTolerance;
        
///////////SET UP GPUs
        set_local_rank(data, &localRank);

        int number_devices = acc_get_num_devices(acc_device_nvidia);

        

        acc_set_device_num(localRank % number_devices,acc_device_nvidia);
        acc_init(acc_device_nvidia);
        int id = acc_get_device_num(acc_device_nvidia);
        
        printf("LOCAL RANK %d DEVICES %d  ASSIGNMENT %d ID-GPU%d \n",localRank, number_devices, (localRank%number_devices),id);

        #pragma acc data pcopyin(uold[0:LENGTH],afF[0:LENGTH]) pcopy(afU[0:LENGTH])               
        {

            while (data->iIterCount < data->iIterMax && residual > data->fTolerance) 
            {

                residual = 0.0;
                residualup = 0.0;
                residualdown = 0.0;
                residualwest = 0.0;
                residualeast = 0.0;
                residualinner = 0.0;
                
                temp = uold;
                uold = afU;
                afU = temp;

                //boundaries nord
                j = data->iRowFirst + 1;

                #pragma acc parallel loop present(uold[LENGTH],afU[LENGTH],afF[LENGTH]) private(j, i, fLRes) reduction(+:residualup)
                for (i = data->iColFirst + 1; i <= data->iColLast - 1; i++)
                {
                        fLRes = ( ax * (UOLD(j, i-1) + UOLD(j, i+1)) + ay * (UOLD(j-1, i) + UOLD(j+1, i)) +  b * UOLD(j, i) - F(j, i)) / b;
                        U(j,i) = UOLD(j,i) - data->fRelax * fLRes;
                        residualup += fLRes * fLRes;
                    
                }
                //boundaries south
                j = data->iRowLast - 1;

                #pragma acc parallel loop present(uold[LENGTH],afU[LENGTH],afF[LENGTH]) private(j, i, fLRes) reduction(+:residualdown)
                for (i = data->iColFirst + 1; i <= data->iColLast - 1; i++)
                {    
                        fLRes = ( ax * (UOLD(j, i-1) + UOLD(j, i+1)) + ay * (UOLD(j-1, i) + UOLD(j+1, i)) +  b * UOLD(j, i) - F(j, i)) / b;
                        U(j,i) = UOLD(j,i) - data->fRelax * fLRes;
                        residualdown += fLRes * fLRes;  
                }
                //boundaries west
                i = data->iColFirst + 1;
                #pragma acc parallel loop present(uold[LENGTH],afU[LENGTH],afF[LENGTH]) private(j, i, fLRes) reduction(+:residualwest)
                for (j = data->iRowFirst + 1; j <= data->iRowLast - 1; j++)
                {

                        fLRes = ( ax * (UOLD(j, i-1) + UOLD(j, i+1)) + ay * (UOLD(j-1, i) + UOLD(j+1, i)) +  b * UOLD(j, i) - F(j, i)) / b;
                        U(j,i) = UOLD(j,i) - data->fRelax * fLRes;
                        residualwest += fLRes * fLRes;
                    
                }
                //boundaries east
                i = data->iColLast - 1;
                #pragma acc parallel loop present(uold[LENGTH],afU[LENGTH],afF[LENGTH]) private(j, i, fLRes) reduction(+:residualeast)
                for (j = data->iRowFirst + 1; j <= data->iRowLast - 1; j++)
                {

                        fLRes = ( ax * (UOLD(j, i-1) + UOLD(j, i+1)) + ay * (UOLD(j-1, i) + UOLD(j+1, i)) +  b * UOLD(j, i) - F(j, i)) / b;
                        U(j,i) = UOLD(j,i) - data->fRelax * fLRes;
                        residualeast += fLRes * fLRes;
                    
                }

                //inner
                #pragma acc parallel loop present(uold[LENGTH],afU[LENGTH],afF[LENGTH]) private(j, i, fLRes) async
                for (j = data->iRowFirst + 2; j <= data->iRowLast - 2; j++)
                {
                    #pragma acc loop reduction(+:residualinner)
                    for (i = data->iColFirst + 2; i <= data->iColLast - 2; i++)
                    {
                        fLRes = ( ax * (UOLD(j, i-1) + UOLD(j, i+1))
                                + ay * (UOLD(j-1, i) + UOLD(j+1, i))
                                +  b * UOLD(j, i) - F(j, i)) / b;

                        //UPDATE SOLUTION
                        U(j,i) = UOLD(j,i) - data->fRelax * fLRes;
                        
                        /* accumulate residual error */
                        residualinner += fLRes * fLRes;
                    }
                }

                #pragma acc data present(uold, afU)
                ExchangeJacobiMpiData(data, uold, afU); 
                
                #pragma acc wait
                /* end for */

            
                tmpResd = residualup+ residualdown+ residualwest+ residualeast+ residualinner;

                MPI_Allreduce(&tmpResd, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                /* error check */
                (data->iIterCount)++;
                residual = sqrt(residual) / (data->iCols * data->iRows);

                

                
            }//while

            data->fResidual = residual;
            
        }
        
        free(uold);
    }
    else 
    {
        fprintf(stderr,"Error: can't allocate memory\n");
        Finish(data);
        exit(1);
    }
}


void ExchangeJacobiMpiData(struct JacobiData *  data, double *  uold, double *  afU){
    //printf("MPI %d CAME BY my proc NULL is %d\n\n\n", data->iMyRank, MPI_PROC_NULL);
    

    //assign communicator for topology
    MPI_Comm topocomm = data->Topo_comm;
    //printf("MPI: %d TOPO COMM %lu\n",data->iMyRank,topocomm);



    int iReqCnt = 0;
    int pos=0;

 

    const int iTagMoveDown = 8;
    const int iTagMoveUp = 9;
    const int iTagMoveWest = 10;
    const int iTagMoveEast = 11;

    
    //double * restrict afF = data->afF;

    int north, west, south, east;
    
    MPI_Cart_shift(topocomm, 0, 1, &west, &east);
    MPI_Cart_shift(topocomm, 1, 1, &north, &south);
    
    //printf("MPI: %d has north %d south %d east %d west %d\n",data->iMyRank,north,south,east,west);

    int width = (data->iColLast - data->iColFirst + 1); 
    int height = (data->iRowLast - data->iRowFirst + 1);
    //printf("WIDTH %d HEIGHT %d\n", width, height);

    double *  toW = acc_malloc((size_t)(height-2)*sizeof(double));
    double *  fW  = acc_malloc((size_t)(height-2)*sizeof(double)); 
    double *  toE = acc_malloc((size_t)(height-2)*sizeof(double));
    double *  fE  = acc_malloc((size_t)(height-2)*sizeof(double)); 
    /*
    double *  toN = acc_malloc((size_t)(width-2)*sizeof(double));
    double *  fN  = acc_malloc((size_t)(width-2)*sizeof(double)); 
    double *  toS = acc_malloc((size_t)(width-2)*sizeof(double));
    double *  fS  = acc_malloc((size_t)(width-2)*sizeof(double));
    */
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //MPI_Request *  request = malloc((size_t)8*sizeof(MPI_Request));
        //#pragma enter data copyin(request)
        //MPI_Status *  status = malloc((size_t)8*sizeof(MPI_Status));
        //#pragma enter data copyin(status)

        //MPI_Request   request[8];
        //MPI_Status   status[8];
        //free(request);
        //free(status);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //printf("MPI: %d, BUFFERS %x %x %x %x %x %x %x %x\n",data->iMyRank,toW, fW, toN, fN, toE, fE, toS, fS);

    #pragma acc data pcreate(uold[LENGTH],toW[(height-2)], fW[(height-2)], toE[(height-2)], fE[(height-2)])//,toN[(width-2)], fN[(width-2)],  toS[(width-2)], fS[(width-2)])  
    {
            //#pragma host_data use_device(toW, fW, toN, fN, toE, fE, toS, fS)
            //printf("MPI gpu: %d, BUFFERS %x %x %x %x %x %x %x %x\n",data->iMyRank,toW, fW, toN, fN, toE, fE, toS, fS);

            if(west != MPI_PROC_NULL){

                    // Send a column
                    #pragma acc parallel loop present(uold[:LENGTH], toW[:(height-2)]) async(1)
                    for(pos=0; pos < (height-2); pos++ ){

                        toW[pos]= UOLD((data->iRowFirst+1)+pos, (data->iColFirst + 1));
                        //printf("toW %lf\n",toW[pos]);
                    }
                    /*
                    #pragma acc host_data use_device(toW) 
                    {
                        printf("MPI %d Swest\n",data->iMyRank);
                        MPI_Isend(&toW,(height-2), MPI_DOUBLE, west, iTagMoveWest, topocomm, &request[iReqCnt]);
                        iReqCnt++;
                    }

                    #pragma acc host_data use_device(fW) 
                    {
                        printf("MPI %d west\n",data->iMyRank);
                        MPI_Irecv(&fW,(height-2), MPI_DOUBLE, west, iTagMoveEast, topocomm,&request[iReqCnt]);
                        iReqCnt++;
                    }*/
                    #pragma acc host_data use_device(toW,fW) 
                    {
                        printf("MPI %d Swest\n",data->iMyRank);
                        MPI_Sendrecv( toW,(height-2),MPI_DOUBLE,west,iTagMoveWest,fW,(height-2),MPI_DOUBLE,west,iTagMoveEast,topocomm,MPI_STATUS_IGNORE );
                    }

                    #pragma acc parallel loop present(uold[:LENGTH], fW[:(height-2)])   
                    for( pos=0; pos < (height-2); pos++ ){

                        UOLD((data->iRowFirst+1)+pos, (data->iColFirst)) = fW[pos];
                        //printf("fW %lf\n",fW[pos]);
                    }
                    
                    printf("DONE WEST \n");
            }

            if(east != MPI_PROC_NULL){
                    // Send a column

                    #pragma acc parallel loop present(uold[:LENGTH],toE[:(height-2)])  
                    for(pos=0; pos < (height-2); pos++ ){

                        toE[pos]= UOLD((data->iRowFirst+1) + pos, (data->iColLast-1));
                        //printf("toE %lf\n",toE[pos]);
                    }
                    /*
                    #pragma acc host_data use_device(toE)  
                    {
                        MPI_Isend(&toE,(height-2), MPI_DOUBLE, east, iTagMoveEast, topocomm, &request[iReqCnt]);
                        iReqCnt++;
                        printf("MPI %d Seast\n",data->iMyRank);
                    }

                    #pragma acc host_data use_device(fE)
                    {
                        MPI_Irecv(&fE,(height-2), MPI_DOUBLE, east, iTagMoveWest, topocomm, &request[iReqCnt]);
                        iReqCnt++;
                        printf("MPI %d east\n",data->iMyRank);
                    }
                    */
                    #pragma acc host_data use_device(toE,fE) 
                    {
                        printf("MPI %d Seast\n",data->iMyRank);
                        MPI_Sendrecv( toE,(height-2),MPI_DOUBLE,east,iTagMoveEast,fE,(height-2),MPI_DOUBLE,east,iTagMoveWest,topocomm,MPI_STATUS_IGNORE );
                    }
                    #pragma acc parallel loop present(uold[:LENGTH], fE[:(height-2)])   
                    for(pos=0; pos < (height-2); pos++ ){

                        UOLD((data->iRowFirst+1) + pos, data->iColLast) = fE[pos];
                        //printf("fE %lf\n",fE[pos]);
                    }
                    printf("DONE EAST \n");
            }

            if(north != MPI_PROC_NULL){
                    /*
                    #pragma acc parallel loop present( uold[:LENGTH], toN[:(width-2)] )  
                    for(pos=0; pos < (width-2); pos++ ){

                        toN[pos]= UOLD( (data->iRowFirst + 1), (data->iColFirst + 1)+pos);
                            
                    }
                    
                    #pragma acc host_data use_device(toN)  
                    {
                        printf("MPI %d Snorth\n",data->iMyRank);

                        MPI_Isend(&toN,  (width-2), MPI_DOUBLE, north, iTagMoveUp, topocomm, &request[iReqCnt]);
                        iReqCnt++;
                    }

                    #pragma acc host_data use_device(fN)
                    {
                        printf("MPI %d north\n",data->iMyRank);

                        MPI_Irecv(&fN,  (width-2), MPI_DOUBLE, north, iTagMoveDown, topocomm, &request[iReqCnt]);
                        iReqCnt++;  
                    }*/
                    #pragma acc host_data use_device(uold) 
                    {
                        printf("MPI %d Snorth\n",data->iMyRank);
                        MPI_Sendrecv( &UOLD( (data->iRowFirst + 1), (data->iColFirst + 1)),(width-2),MPI_DOUBLE,north,iTagMoveUp,&UOLD( (data->iRowFirst), (data->iColFirst + 1)),(width-2),MPI_DOUBLE,north,iTagMoveDown,topocomm,MPI_STATUS_IGNORE );
                    }
                    /*
                    #pragma acc parallel loop present(uold[:LENGTH], fN[:(width-2)])   
                    for(pos=0; pos < (width-2); pos++ ){

                        UOLD( (data->iRowFirst), (data->iColFirst + 1)+ pos) = fN[pos];

                    }*/

            }
            if(south != MPI_PROC_NULL){
                    /*
                    #pragma acc parallel loop present(uold[:LENGTH], toS[:(width-2)])   
                    for(pos=0; pos < (width-2); pos++ ){

                        toS[pos]= UOLD( (data->iRowLast - 1), (data->iColFirst+1)+pos);
                            
                    }
                    
                    #pragma acc host_data use_device(toS)  
                    {
                        printf("MPI %d Ssouth\n",data->iMyRank);

                        MPI_Isend(&toS, (width-2), MPI_DOUBLE, south, iTagMoveDown, topocomm, &request[iReqCnt]);
                        iReqCnt++;
                    }
                    
                    #pragma acc host_data use_device(fS)
                    {
                        printf("MPI %d south\n",data->iMyRank);

                        MPI_Irecv(&fS, (width-2), MPI_DOUBLE, south, iTagMoveUp, topocomm, &request[iReqCnt]);
                        iReqCnt++;
                    }
                    */
                    #pragma acc host_data use_device(uold) 
                    {
                        printf("MPI %d Ssouth\n",data->iMyRank);
                        MPI_Sendrecv( &UOLD( (data->iRowLast - 1), (data->iColFirst+1)),(width-2),MPI_DOUBLE,south,iTagMoveDown,&UOLD( (data->iRowLast), (data->iColFirst+1)),(width-2),MPI_DOUBLE,south,iTagMoveUp,topocomm,MPI_STATUS_IGNORE );
                    }
                    /*
                    #pragma acc parallel loop present(uold[:LENGTH], fS[:(width-2)])  
                    for(pos=0; pos < (width-2); pos++ ){

                        UOLD( (data->iRowLast), (data->iColFirst+1)+pos) = fS[pos];

                    }*/	
            }
    //#pragma acc update host(request)
    //#pragma acc update host(status)  
    }
   
    printf("WAITALL %d \n", iReqCnt);

    //MPI_Waitall(iReqCnt, request, MPI_STATUSES_IGNORE);



    printf("PAST WAITALL\n");
    //#pragma acc exit data delete(request, status)
    acc_free(toE);
    acc_free(fE);
    acc_free(toW);
    acc_free(fW);
    /*
    acc_free(toN);
    acc_free(fN);
    acc_free(toS);
    acc_free(fS);
    */


    printf("MPI %d EXITS FROM EXCHANGE\n",data->iMyRank);
}



void set_local_rank(struct JacobiData *data, int * local){
    int rank = data->iMyRank;
    
    MPI_Comm local_comm;
    MPI_Info info;
    MPI_Info_create(&info);

    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, info, &local_comm);

    
    MPI_Comm_rank(local_comm, local);
    MPI_Comm_free(&local_comm);
    MPI_Info_free(&info);
}

