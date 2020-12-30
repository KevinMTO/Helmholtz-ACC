/* 
**************************************************************************
*                                                                        * 
* Program to solve a finite difference discretization of the             *
* Helmholtz equation:   (d2/dx2)u + (d2/dy2)u - alpha u = f              *
* using Jacobi iterative method, with Dirichlect boundary conditions.    *
*                                                                        * 
* Used as an exercise in the course:                                     *
*                     Introduction to Hybrid Programming in HPC          *
*                                                                        * 
* Exercise = the MPI-only version of the code                            *
*   TODO   = add OpenMP directives to get a hybrid MPI+OpenMP code !     *
*          + do various improvements to achieve better performance !     *
*                                                                        * 
* Provided = Exercise (MPI-only) + various Solutions (different levels)  *


#include <mpi.h>
#include<omp.h>

#if defined(WIN32) || defined (_CLUSTER_OPENMP)
/* fix a visualstudio 2005 bug */
#include <omp.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "realtime.h"
#include "jacobi.h"

#define U(j,i) data->afU[((j) - data->iRowFirst) * data->iCols + (i)]
#define F(j,i) data->afF[((j) - data->iRowFirst) * data->iCols + (i)]

/*
 * setting values, init mpi, omp etc
 */
void Init(struct JacobiData *data, int *argc, char **argv)
{
    int i;
    int block_lengths[8];
    int provided;
    MPI_Datatype MY_JacobiData;
    MPI_Datatype typelist[8] = { MPI_INT,  MPI_INT, MPI_INT, MPI_INT, MPI_INT, 
          MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };

    MPI_Aint displacements[8];

    /* MPI Initialization */
    //############################################################################
    //MPI_Init(argc, &argv); //only MPI
    MPI_Init_thread(argc, &argv, MPI_THREAD_FUNNELED, &provided); //MPI+ OPENMP
    //#########################################################################
    MPI_Comm_rank(MPI_COMM_WORLD, &data->iMyRank);
    MPI_Comm_size(MPI_COMM_WORLD, &data->iNumProcs);

    if (data->iMyRank == 0)
    {
/* default medium */
        data->iCols      = 2000;
        data->iRows      = 2000;
        data->fAlpha     = 0.8;
        data->fRelax     = 1.0;
        data->fTolerance = 1e-10;
        data->iIterMax   = 50;
#ifdef READ_INPUT
        printf("Input n - matrix size in x direction:                 ");
        scanf("%d", &data->iCols);
        printf("\nInput m - matrix size in y direction:               ");
        scanf("%d", &data->iRows);
        printf("\nInput alpha - Helmholtz constant:                   ");
        scanf("%lf", &data->fAlpha);
        printf("\nInput relax - Successive over-relaxation parameter: ");
        scanf("%lf", &data->fRelax);
        printf("\nInput tol - error tolerance for iterrative solver:  ");
        scanf("%lf", &data->fTolerance);
        printf("\nInput mits - Maximum iterations for solver:         ");
        scanf("%d", &data->iIterMax);
#elif defined DATA_LARGE
        data->iCols      = 7000;
        data->iRows      = 7000;
        data->fAlpha     = 0.8;
        data->fRelax     = 1.0;
        data->fTolerance = 1e-12;
        data->iIterMax   = 2;
#elif defined DATA_SMALL
        data->iCols      = 200;
        data->iRows      = 200;
        data->fAlpha     = 0.8;
        data->fRelax     = 1.0;
        data->fTolerance = 1e-7;
        data->iIterMax   = 1000;
#endif
        printf("\n-> matrix size: %dx%d"
               "\n-> alpha: %f"
               "\n-> relax: %f"
               "\n-> tolerance: %f"
               "\n-> #of iterations: %d \n\n",
               data->iCols, data->iRows, data->fAlpha, data->fRelax,
               data->fTolerance, data->iIterMax);
    }

    /* Build MPI Datastructure */
    for(i = 0; i < 8; i++)
    {
         block_lengths[i] = 1;
    }
    displacements [0] = (MPI_Aint)offsetof(struct JacobiData, iRows);
    displacements [1] = (MPI_Aint)offsetof(struct JacobiData, iCols);
    displacements [2] = (MPI_Aint)offsetof(struct JacobiData, iRowFirst);
    displacements [3] = (MPI_Aint)offsetof(struct JacobiData, iRowLast);
    displacements [4] = (MPI_Aint)offsetof(struct JacobiData, iIterMax);
    displacements [5] = (MPI_Aint)offsetof(struct JacobiData, fAlpha);
    displacements [6] = (MPI_Aint)offsetof(struct JacobiData, fRelax);
    displacements [7] = (MPI_Aint)offsetof(struct JacobiData, fTolerance);
    
    MPI_Type_create_struct(8, block_lengths, displacements, typelist, &MY_JacobiData);
    MPI_Type_commit(&MY_JacobiData);
    
    /* Send input parameters to all procs */
    MPI_Bcast(data, 1, MY_JacobiData, 0, MPI_COMM_WORLD);

    /* calculate bounds for the task working area */
    data->iRowFirst = data->iMyRank * (data->iRows - 2) / data->iNumProcs;
    if (data->iMyRank == data->iNumProcs - 1) 
    {
        data->iRowLast = data->iRows - 1;
    }
    else
    {
        data->iRowLast =
            (data->iMyRank + 1) * (data->iRows - 2) / data->iNumProcs + 1;
        //printf("MPI: %d LAST: %d--myr %d, ir: %d , inum+1 : %d\n\n",data->iMyRank,data->iRowLast,(data->iMyRank + 1), (data->iRows - 2), (data->iNumProcs + 1));
        
    }
    
    data->afU = (double*) malloc(
        (data->iRowLast - data->iRowFirst + 1) * data->iCols * sizeof(double));
    data->afF = (double*) malloc(
        (data->iRowLast - data->iRowFirst + 1) * data->iCols * sizeof(double));

    /* calculate dx and dy */
    data->fDx = 2.0 / (data->iCols - 1);
    data->fDy = 2.0 / (data->iRows - 1);

    data->iIterCount = 0;
    //printf("MPI: %d has first: %d and last %d on total : %d \n\n",data->iMyRank,data->iRowFirst,data->iRowLast, data->iRows);
    return;
}

/*
 * final cleanup routines
 */
void Finish(struct JacobiData * data)
{
    free (data->afU);
    free (data->afF);
    MPI_Finalize();

    return;
}

/*
 * print result summary
 */
void PrintResults(const struct JacobiData * data)
{
    if (data->iMyRank == 0)
    {
        printf(" Number of iterations : %d\n", data->iIterCount);
        printf(" Residual             : %le\n", data->fResidual);
        printf(" Solution Error       : %1.12lf\n", data->fError);
        printf(" Elapsed Time         : %5.7lf\n", 
               data->fTimeStop-data->fTimeStart);
        printf(" MFlops               : %6.6lf\n", 
            0.000013 * data->iIterCount * (data->iCols - 2) * (data->iRows - 2)
            / (data->fTimeStop - data->fTimeStart));
    }

    return;
}

/*
 * Initializes matrix
 * Assumes exact solution is u(x,y) = (1-x^2)*(1-y^2)
 */
void InitializeMatrix(struct JacobiData * data)
{
    int i, j;
    double xx, yy, xx2, yy2;

    /* Initialize initial condition and RHS */
    #pragma omp parallel for private(i, j,xx, yy, xx2, yy2)
    for (j = data->iRowFirst; j <= data->iRowLast; j++)
    {
        for (i = 0; i < data->iCols; i++)
        {
            /* TODO: check if this values have to be ints or doubles */
            xx = (double) (-1.0 + data->fDx * i);
            yy = (double) (-1.0 + data->fDy * j);

            xx2 = xx * xx;
            yy2 = yy * yy;

            U(j,i) = 0.0;
            F(j,i) = -data->fAlpha * (1.0 - xx2) * (1.0 - yy2)
                     + 2.0 * (-2.0 + xx2 + yy2);

            //if(i==0){
            //printf("th: %d, i: %d, j: %d\n", omp_get_thread_num(),i,j);
            //}
        }
    }
}

/*
 * Checks error between numerical and exact solution
 */
void CheckError(struct JacobiData * data)
{
    double error = 0.0;
    int i, j;
    double xx, yy, temp;
    double local_error= 0.0;

    #pragma omp parallel for private(i, j, xx, yy, local_error, temp)
    for (j = data->iRowFirst; j <= data->iRowLast; j++)
    {
        if ((data->iMyRank != 0 && j == data->iRowFirst) || 
            (data->iMyRank != data->iNumProcs - 1 && j == data->iRowLast))
            continue;

        for (i = 0; i < data->iCols; i++)
        {
            xx   = -1.0 + data->fDx * i;
            yy   = -1.0 + data->fDy * j;
            temp = U(j,i) - (1.0 - xx * xx) * (1.0 - yy * yy);
            local_error += temp * temp;

            //if(i==0){
            //printf("th: %d, i: %d, j: %d\n", omp_get_thread_num(),i,j);
            //}
        }
        #pragma omp critical
		{
			error+=local_error;
    	}
    }

    data->fError = error;
    MPI_Reduce(&data->fError, &error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    data->fError = sqrt(error) / (data->iCols * data->iRows);
        
    return;
}

int main (int argc, char** argv)
{
    int retVal = 0;    /* return value */

    struct JacobiData myData;

    /* sets default values or reads from stdin
     * inits MPI and OpenMP if needed
     * distribute MPI data, calculate MPI bounds
     */
    Init(&myData, &argc, argv);

    if (myData.afU && myData.afF)
    {
        /* matrix init */
        InitializeMatrix(&myData);

        /* starting timer */
        myData.fTimeStart = GetRealTime();

        /* running calculations */
        Jacobi(&myData);

        /* stopping timer */
        myData.fTimeStop = GetRealTime();

        /* error checking */
        CheckError(&myData);

        /* print result summary */
        PrintResults(&myData);
    }
    else
    {
        printf(" Memory allocation failed ...\n");
        retVal = -1;
    }

    /* cleanup */
    Finish(&myData);
    
    return retVal;
}

