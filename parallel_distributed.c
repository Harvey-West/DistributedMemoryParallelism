#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* A two dimensional array can be expressed by one dimensional malloc whose size
is squared.
To enable two dimensional access, simply malloc pointers to each
row index.
*/
double **createArrayOfDoubles(int numberOfRows, int sizeOfRow);
/* Fills each boundary index with 1.0
*/
void populateBoundaryValues(double **array, int sizeOfRow, int startYIndex,
                            int endYIndex, int numberOfVisibleRows);
/*Returns the number of rows that a thread will work on, i.e. the grain.*/
int getGrainSize(int sizeOfRow, int numberOfThreads, int myRank);
/*Although a thread will work on the grain, it will need to see the boundary 
values. This may mean seeing 1 additional row or 2. This function returns
how many rows each thread will need to use to work*/
int getNumberOfVisibleRows(int myRank, int grainSize, int nproc);
void printArray(double **array, int sizeOfRow);
void printGrain(double **array, int offset, int grainSize, int sizeOfRow,
                int nproc);
void printVisibleArea(double **array, int numberOfVisibleRows, int sizeOfRow);

int main(int argc, char **argv)
{
    int sizeOfRow = strtol(argv[1], NULL, 10),
        numberOfThreads = strtol(argv[2], NULL, 10),
        grainSize,
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

    MPI_Status stat;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    grainSize = getGrainSize(sizeOfRow, numberOfThreads, myRank);
    int startYIndex = myRank * grainSize,
        endYIndex = startYIndex + (grainSize - 1);
    double **workingArray, **destinationArray, **wholeArray;

    printf("hello world from ’%d’, grainSize = %d\n", myRank, grainSize);

    int numberOfVisibleRows = getNumberOfVisibleRows(myRank, grainSize, nproc);

    workingArray = createArrayOfDoubles(numberOfVisibleRows, sizeOfRow);
    destinationArray = createArrayOfDoubles(numberOfVisibleRows, sizeOfRow);

    populateBoundaryValues(workingArray, sizeOfRow, startYIndex, endYIndex,
                           numberOfVisibleRows);
    populateBoundaryValues(destinationArray, sizeOfRow, startYIndex, endYIndex,
                           numberOfVisibleRows);

    if (myRank == 0)
    {
        wholeArray = createArrayOfDoubles(sizeOfRow, sizeOfRow);
        populateBoundaryValues(wholeArray, sizeOfRow, startYIndex,
                               sizeOfRow - 1, sizeOfRow);
    }

    //TODO
    /*
    1)  Get one pass then send back to 0 working.
    2)  Get one pass then broadcast if not relaxed, get all other threads
        to stop printing some statement once this is done.
    2)  Each thread performs calculation, while calculating check if not relaxed.
        a)  If not relaxed, broadcast it is not relaxed to all other threads.
            All threads stop checking on that number of passes if relaxed.
    3)  Finished thread calculation and result was it was relaxed.
        b)  If all say relaxed and finished, reduce.
        c)  If relaxed and nearest neighbours 
    if finish, broadcast relaxed. Every thread keeps count of number of passes.
    When broadcasting if relaxed, make it a tuple of isRelaxed and number of 
    passes.
    */
    //------------------------

    if (myRank == 0)
    {
        int x;
        for (x = 0; x < sizeOfRow * grainSize; x++)
        {
            *(wholeArray[1] + x) = *(workingArray[1] + x);
        }
        /*Noting that MPI_Recv asks for a number of items to receive.
        Importantly it specifies the *maximum* number of items. As such each, if
        carefully programmed so that no other send uses the same tag I can 
        use a larger than actual maximum count here to encapsulate every
        receive from every thread. No data will overlap as each send is careful
        with what data it sends and how much.
        */
        int threadIndex,
            maximumAmount = (grainSize+1) * sizeOfRow;
        for (threadIndex = 1; threadIndex < nproc - 1; threadIndex++)
        {
            MPI_Recv(wholeArray[threadIndex * grainSize], maximumAmount,
                     MPI_DOUBLE, threadIndex, threadIndex, MPI_COMM_WORLD,
                     &stat);
        }
    }
    else
    {
        MPI_Send(workingArray[1], grainSize * sizeOfRow,
                 MPI_DOUBLE, 0, myRank, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* implicit barrier in Finalize */
    /*MPI_Barrier(MPI_COMM_WORLD);*/
    MPI_Finalize();
    if (myRank == 0)
    {
        printArray(wholeArray, sizeOfRow);
        // printGrain(workingArray, myRank, grainSize, sizeOfRow, nproc);
        // printVisibleArea(workingArray, numberOfVisibleRows, sizeOfRow);
    }
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
        for (x = 0; x < sizeOfRow - 1; x++)
        {
            printf("%f,", array[y][x]);
        }
        printf("%f\n", array[y][x]);
    }
    printf("\n");
}
void printGrain(double **array, int myRank, int grainSize, int sizeOfRow,
                int nproc)
{
    //Would need a barrier to prevent other threads printing.
    int y = 1, endYIndex, x;
    if (myRank == 0 || myRank == nproc - 1)
    {
        endYIndex = grainSize;
    }
    else
    {
        endYIndex = grainSize + 1;
    }
    for (y; y < endYIndex; y++)
    {
        printf("Y = %d  ,", myRank * grainSize + y);
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
int getGrainSize(int sizeOfRow, int numberOfThreads, int myRank)
{
    int grainSize = 0, numberOfViableRows = sizeOfRow - 2,
        remaining = numberOfViableRows % numberOfThreads;
    if (remaining)
    {
        grainSize = (numberOfViableRows - remaining) / numberOfThreads;
        //Spreads work across threads, aside from 0 as 0 does more other work.
        if ((remaining - myRank) > -1 && myRank)
        {
            grainSize += 1;
        }
    }
    else
    {
        grainSize = numberOfViableRows / numberOfThreads;
    }
    return grainSize;
}