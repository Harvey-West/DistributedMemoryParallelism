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
    //Initialise all values to 0
    double *array = calloc(sizeOfArray * sizeOfArray, sizeof(double *));
    double **rowPointers = malloc(sizeOfArray * sizeof(double *));
    for (rowIndex = 0; rowIndex < sizeOfArray; rowIndex++)
    {
        rowPointers[rowIndex] = &array[rowIndex * sizeOfArray];
    }
    return rowPointers;
}
/* Fills each index with random double from 0.0 to 1 (not inclusive).
*/
void populateBoundaryValues(double **array, int sizeOfArray, int startYIndex,
                            int endYIndex)
{
    int y, x, endXIndex = sizeOfArray - 1;

    for (y = startYIndex; y <= endYIndex; y++)
    {
        array[y][0] = 1.0;
        array[y][endXIndex] = 1.0;
    }
    if (startYIndex == 0)
    {
        for (x = 1; x < endXIndex; x++)
        {
            array[0][x] = 1.0;
        }
    }
    if (endYIndex == endXIndex)
    {
        for (x = 1; x < endXIndex; x++)
        {
            array[endYIndex][x] = 1.0;
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
    int sizeOfRow = strtol(argv[1], NULL, 10),
        numberOfThreads = strtol(argv[2], NULL, 10),
        countOfPasses = 0;
    double userPrecision;
    sscanf(argv[3], "%lf", &userPrecision);
    int rc, myrank, nproc, namelen;
    char name[MPI_MAX_PROCESSOR_NAME];
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
        double **array = createArrayOfDoubles(sizeOfRow);
        populateBoundaryValues(array, sizeOfRow, 0, sizeOfRow-1);
        printArray(array, sizeOfRow);
    }
    namelen = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(name, &namelen);
    printf("hello world %d from ’%s’\n", myrank, name);
    /* implicit barrier in Finalize */
    /*MPI_Barrier(MPI_COMM_WORLD);*/
    MPI_Finalize();
    return 0;
}