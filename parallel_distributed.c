#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* A two dimensional array can be expressed by one dimensional malloc whose size
is squared.
To enable two dimensional access, simply malloc pointers to each
row index.
*/
double **createArrayOfDoubles(int numberOfRows, int sizeOfRow)
{
    int rowIndex;
    //Initialise all values to 0
    double *array = calloc(numberOfRows * sizeOfRow, sizeof(double *));
    double **rowPointers = malloc(sizeOfRow * sizeof(double *));
    for (rowIndex = 0; rowIndex < numberOfRows; rowIndex++)
    {
        rowPointers[rowIndex] = &array[rowIndex * sizeOfRow];
    }
    return rowPointers;
}
/* Fills each index with random double from 0.0 to 1 (not inclusive).
*/
void populateBoundaryValues(double **array, int sizeOfRow, int startYIndex,
                            int endYIndex, int numberOfRowsVisible)
{
    int y, x, endXIndex = sizeOfRow - 1;

    for (y = 0; y < numberOfRowsVisible; y++)
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
            array[numberOfRowsVisible - 1][x] = 1.0;
        }
    }
}
void printArray(double **array, int sizeOfRow);
void printGrain(double **array, int offset, int grainSize, int sizeOfRow);
void printVisibleArea(double **array, int numberOfRowsVisible, int sizeOfRow);

int main(int argc, char **argv)
{
    int sizeOfRow = strtol(argv[1], NULL, 10),
        numberOfThreads = strtol(argv[2], NULL, 10),
        grainSize = sizeOfRow / numberOfThreads,
        countOfPasses = 0;
    double userPrecision;
    sscanf(argv[3], "%lf", &userPrecision);
    int rc, myRank, nproc, namelen;
    char name[MPI_MAX_PROCESSOR_NAME];
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        printf("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    int startYIndex = myRank * grainSize,
        endYIndex = startYIndex + (grainSize - 1);
    double **array;
    printf("hello world from ’%d’\n", myRank);

    int numberOfRowsVisible = 0;
    if (myRank == 0 || myRank == nproc - 1)
    {
        numberOfRowsVisible = grainSize + 1;
    }
    else
    {
        numberOfRowsVisible = grainSize + 2;
    }
    array = createArrayOfDoubles(numberOfRowsVisible, sizeOfRow);
    populateBoundaryValues(array, sizeOfRow, startYIndex, endYIndex,
                           numberOfRowsVisible);

    printGrain(array, myRank, grainSize, sizeOfRow);
    // printVisibleArea(array,numberOfRowsVisible,sizeOfRow);
    
    /* implicit barrier in Finalize */
    /*MPI_Barrier(MPI_COMM_WORLD);*/
    MPI_Finalize();
    return 0;
}


void printArray(double **array, int sizeOfRow)
{
    int y, x;
    for (y = 0; y < sizeOfRow; y++)
    {
        for (x = 0; x < sizeOfRow; x++)
        {
            printf("%f,", array[y][x]);
        }
        printf("\n");
    }
    printf("\n");
}
void printGrain(double **array, int offset, int grainSize, int sizeOfRow)
{
    //Would need a barrier to prevent other threads printing.
    int y = 0, endYIndex, x;
    if (offset == 0)
    {
        endYIndex = grainSize;
    }
    else
    {
        y = 1;
        endYIndex = grainSize + 1;
    }
    for (y; y < endYIndex; y++)
    {
        for (x = 0; x < sizeOfRow; x++)
        {
            printf("%f,", array[y][x]);
        }
        printf("\n");
    }
    printf("\n");
}
void printVisibleArea(double **array, int numberOfRowsVisible, int sizeOfRow)
{
    //Would need a barrier to prevent other threads printing.
    int y, x;
    for (y = 0; y < numberOfRowsVisible; y++)
    {
        for (x = 0; x < sizeOfRow; x++)
        {
            printf("%f,", array[y][x]);
        }
        printf("\n");
    }
    printf("\n");
}