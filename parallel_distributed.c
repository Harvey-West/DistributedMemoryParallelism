#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* A two dimensional array can be expressed by one dimensional malloc whose size
is squared.
To enable two dimensional access, simply malloc pointers to each
row index.
*/
double **createArrayOfDoubles(int sizeOfArray)
{
    int rowIndex;
    double *array = malloc(sizeOfArray * sizeOfArray * sizeof(double *));
    double **rowPointers = malloc(sizeOfArray * sizeof(double*));
    for (rowIndex = 0; rowIndex < sizeOfArray; rowIndex++)
    {
        rowPointers[rowIndex] = &array[rowIndex*sizeOfArray];
    }
    return rowPointers;
}
/* Fills each index with random double from 0.0 to 1 (not inclusive).
*/
void populateBoundaryValues(double **array, int sizeOfArray)
{
    int y, x;

    for (y = 0; y < sizeOfArray; y++)
    {
        for (x = 0; x < sizeOfArray; x++)
        {
            double randNum = ((double)rand() / (double)RAND_MAX);
            array[y][x] = randNum;
        }
    }
}
void printArray(double **array, int sizeOfArray)
{
    int y, x;
    for (y = 0; y < sizeOfArray; y++)
    {
        for (x = 0; x < sizeOfArray; x++)
        {
            printf("%f,", array[y][x]);
        }
        printf("\n");
    }
    printf("\n");
}
int main(int argc, char **argv)
{
    int rc, myrank, nproc, namelen;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Status stat;
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        printf("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if (myrank == 0)
    {
        MPI_Recv(n, 1, MPI_INT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &stat);
        printf("Hello from 1 to 0: %d\n", *n);
    }
    namelen = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(name, &namelen);
    printf("hello world %d from ’%s’\n", myrank, name);
    /* implicit barrier in Finalize */
    /*MPI_Barrier(MPI_COMM_WORLD);*/
    MPI_Finalize();
    return 0;
}