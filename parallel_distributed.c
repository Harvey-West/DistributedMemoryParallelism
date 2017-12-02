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
                            int endYIndex, int numberOfVisibleRows)
{
    int y, x, endXIndex = sizeOfRow - 1;

    for (y = 0; y < numberOfVisibleRows; y++)
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
            array[numberOfVisibleRows - 1][x] = 1.0;
        }
    }
}
int getNumberOfVisibleRows(int myRank, int grainSize, int nproc);
void printArray(double **array, int sizeOfRow);
void printGrain(double **array, int offset, int grainSize, int sizeOfRow);
void printVisibleArea(double **array, int numberOfVisibleRows, int sizeOfRow);

int main(int argc, char **argv)
{
    int sizeOfRow = strtol(argv[1], NULL, 10),
        numberOfThreads = strtol(argv[2], NULL, 10),
        grainSize = sizeOfRow / numberOfThreads,
        countOfPasses = 0;
    double userPrecision;
    sscanf(argv[3], "%lf", &userPrecision);
    int rc, myRank, nproc;

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
    double **workingArray, **destinationArray, **wholeArray;

    printf("hello world from ’%d’\n", myRank);

    int numberOfVisibleRows = getNumberOfVisibleRows(myRank, grainSize, nproc);

    workingArray = createArrayOfDoubles(numberOfVisibleRows, sizeOfRow);
    destinationArray = createArrayOfDoubles(numberOfVisibleRows, sizeOfRow);

    populateBoundaryValues(workingArray, sizeOfRow, startYIndex, endYIndex,
                           numberOfVisibleRows);

    // printGrain(array, myRank, grainSize, sizeOfRow);
    // printVisibleArea(destinationArray,numberOfVisibleRows,sizeOfRow);

    if (myRank == 0)
    {
        wholeArray = createArrayOfDoubles(sizeOfRow, sizeOfRow);
        populateBoundaryValues(wholeArray, sizeOfRow, 0, sizeOfRow - 1,
                               sizeOfRow);
        printVisibleArea(wholeArray, sizeOfRow, sizeOfRow);
    }
    /* implicit barrier in Finalize */
    /*MPI_Barrier(MPI_COMM_WORLD);*/
    MPI_Finalize();
    return 0;
}

int getNumberOfVisibleRows(int myRank, int grainSize, int nproc)
{
    int numberOfVisibleRows = 0;
    if (myRank == 0 || myRank == nproc - 1)
    {
        numberOfVisibleRows = grainSize + 1;
    }
    else
    {
        numberOfVisibleRows = grainSize + 2;
    }
    return numberOfVisibleRows;
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
void printVisibleArea(double **array, int numberOfVisibleRows, int sizeOfRow)
{
    //Would need a barrier to prevent other threads printing.
    int y, x;
    for (y = 0; y < numberOfVisibleRows; y++)
    {
        for (x = 0; x < sizeOfRow; x++)
        {
            printf("%f,", array[y][x]);
        }
        printf("\n");
    }
    printf("\n");
}