

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
#include "jacobi_enhance.h"

#define U(j,i) data->afU[((j) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + (i- data->iColFirst)]
#define F(j,i) data->afF[((j) - data->iRowFirst) * (data->iColLast - data->iColFirst+1) + (i- data->iColFirst)]





/*
 * setting values, init mpi, omp etc
 */
void Init(struct JacobiData *data, int *argc, char **argv)
{

    int provided;


    /* MPI Initialization */
    //############################################################################
    //MPI_Init(argc, &argv); //only MPI
    MPI_Init_thread(argc, &argv, MPI_THREAD_FUNNELED, &provided); //MPI+ OPENMP
    //#########################################################################
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
        printf("Input n - matrix size in x direction:  ");
        scanf("%d", &data->iCols);
        printf("\nInput m - matrix size in y direction: ");
        scanf("%d", &data->iRows);
        printf("\nInput alpha - Helmholtz constant:  ");
        scanf("%lf", &data->fAlpha);
        printf("\nInput relax - Successive over-relaxation parameter: ");
        scanf("%lf", &data->fRelax);
        printf("\nInput tol - error tolerance for iterrative solver:  ");
        scanf("%lf", &data->fTolerance);
        printf("\nInput mits - Maximum iterations for solver:  ");
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
	//WHAT I'M WORKING WITH
        printf("\n-> matrix size: %dx%d"
               "\n-> alpha: %f"
               "\n-> relax: %f"
               "\n-> tolerance: %f"
               "\n-> #of iterations: %d \n\n",
               data->iCols, data->iRows, data->fAlpha, data->fRelax,
               data->fTolerance, data->iIterMax);
    }


    /* Build MPI Datastructure */
    int i;
    MPI_Datatype MY_JacobiData;

    int block_lengths[10];
    MPI_Datatype typelist[10] = { MPI_INT,  MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint displacements[10];

    for(i = 0; i < 10; i++){
         block_lengths[i] = 1;
    }

    displacements [0] = (MPI_Aint)offsetof(struct JacobiData, iRows);
    displacements [1] = (MPI_Aint)offsetof(struct JacobiData, iCols);
    displacements [2] = (MPI_Aint)offsetof(struct JacobiData, iRowFirst);
    displacements [3] = (MPI_Aint)offsetof(struct JacobiData, iRowLast);
    displacements [4] = (MPI_Aint)offsetof(struct JacobiData, iColFirst);
    displacements [5] = (MPI_Aint)offsetof(struct JacobiData, iColLast);
    displacements [6] = (MPI_Aint)offsetof(struct JacobiData, iIterMax);
    displacements [7] = (MPI_Aint)offsetof(struct JacobiData, fAlpha);
    displacements [8] = (MPI_Aint)offsetof(struct JacobiData, fRelax);
    displacements [9] = (MPI_Aint)offsetof(struct JacobiData, fTolerance);


    
    MPI_Type_create_struct(10, block_lengths, displacements, typelist, &MY_JacobiData);
    MPI_Type_commit(&MY_JacobiData);
    
    /* Send input parameters to all procs */
    MPI_Bcast(data, 1, MY_JacobiData, 0, MPI_COMM_WORLD);



///////////////////////////////////////////////////////////////////////////////////////////////////////////
/*CREATE 2D TOPOLOGY*/
    
    int dimensions[2]={0, 0};

    MPI_Dims_create(data->iNumProcs, 2, dimensions);
    int proc_X = dimensions[0];
    int proc_Y = dimensions[1];

    int periods[2] = {0,0}; //no periods
    

    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, 0, &(data->Topo_comm)); // no reorder for now let's see

    int coords[2];
    MPI_Cart_coords(data->Topo_comm, data->iMyRank, 2, coords);

    int offsetx = coords[0];
    int offsety = coords[1];

///////////////////////////////////////////////////////////////////////////////////////////////////////////

    /* calculate bounds for the task working area */


    data->iRowFirst = offsety * (data->iRows - 2) / proc_Y;
    data->iColFirst = offsetx * (data->iCols - 2) / proc_X;
   

    if (offsety == proc_Y - 1) 
    {
        data->iRowLast = data->iRows - 1;
    }
    else
    {
        data->iRowLast =
            (offsety + 1) * (data->iRows - 2) / proc_Y + 1;;

    }
    if (offsetx == proc_X - 1) 
    {
        data->iColLast = data->iCols - 1;
    }
    else
    {
        data->iColLast =
            (offsetx + 1) * (data->iCols - 2) / proc_X + 1;

    }
    printf("MPI: %d RF: %d RL %d CF %d CL %d -- ofx: %d ofy: %d ir: %d , ic : %d \n\n",data->iMyRank, data->iRowFirst, data->iRowLast, data->iColFirst, data->iColLast, 
     offsetx, offsety, data->iRowLast - data->iRowFirst, data->iColLast - data->iColFirst  );   
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    data->afU = (double*) malloc(
        (data->iRowLast - data->iRowFirst + 1) * (data->iColLast - data->iColFirst+1) * sizeof(double));
    data->afF = (double*) malloc(
        (data->iRowLast - data->iRowFirst + 1) * (data->iColLast - data->iColFirst+1) * sizeof(double));

    /* calculate dx and dy */
    data->fDx = 2.0 / (data->iCols - 1);
    data->fDy = 2.0 / (data->iRows - 1);

    data->iIterCount = 0;

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
        for (i = data->iColFirst; i <= data->iColLast; i++)
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

        for (i = data->iColFirst; i <= data->iColLast; i++)
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
        printf("FINISHED \n\n\n");
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

