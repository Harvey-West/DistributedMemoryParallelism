#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
struct RelaxationVariables
{
    int sizeOfRow, grainSize;
    double precision;
    double **originalArrayPointer, **destinationArrayPointer;
};
/* A two dimensional array can be expressed by one dimensional malloc whose size
is squared.
To enable two dimensional access, simply malloc pointers to each
row index.
*/
double **createArrayOfDoubles(int numberOfRows, int sizeOfRow);
/* Fills each boundary index with 1.0
*/
void populateBoundaryValues(double **array, int sizeOfRow, int startYIndex,
                            int lastVisibleYIndex, int numberOfVisibleRows);
/*Returns the number of rows that a thread will work on, i.e. the grain.*/
int getGrainSize(int sizeOfRow, int numberOfThreads, int myRank);
/*Although a thread will work on the grain, it will need to see the boundary 
values. This may mean seeing 1 additional row or 2. This function returns
how many rows each thread will need to use to work*/
int getNumberOfVisibleRows(int myRank, int grainSize, int nproc);

void performRelaxation(struct RelaxationVariables);

void printArray(double **array, int sizeOfRow);
void printGrain(double **array, int offset, int grainSize, int sizeOfRow,
                int nproc);
void printVisibleArea(double **array, int numberOfVisibleRows, int sizeOfRow);

int main(int argc, char **argv)
{
    int sizeOfRow = strtol(argv[1], NULL, 10),
        numberOfThreads = strtol(argv[2], NULL, 10),
        grainSize,
        remainder = (sizeOfRow - 2) % numberOfThreads,
        countOfPasses = 0,
        rc, myRank, nproc;
    double userPrecision;
    sscanf(argv[3], "%lf", &userPrecision);

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
    int firstVisibleYIndex,
        numberOfVisibleRows = grainSize + 2,
        lastVisibleYIndex;
    if (remainder)
    {
        if (myRank > remainder)
        {
            firstVisibleYIndex = myRank * grainSize + remainder;
            lastVisibleYIndex = firstVisibleYIndex + grainSize + 1;
        }
        else if (myRank == 0)
        {
            firstVisibleYIndex = myRank * grainSize;
            lastVisibleYIndex = firstVisibleYIndex + grainSize + 1;
        }
        else
        {
            firstVisibleYIndex = myRank * grainSize - 1;
            lastVisibleYIndex = firstVisibleYIndex + grainSize + 1;
        }
    }
    else
    {
        firstVisibleYIndex = myRank * grainSize;
        lastVisibleYIndex = firstVisibleYIndex + grainSize + 1;
    }

    double **workingArray, **destinationArray, **wholeArray, **tempPointer;

    printf("hello world from ’%d’, grainSize = %d, startY = %d\n", myRank, grainSize, firstVisibleYIndex);
    MPI_Barrier(MPI_COMM_WORLD);

    workingArray = createArrayOfDoubles(numberOfVisibleRows, sizeOfRow);
    destinationArray = createArrayOfDoubles(numberOfVisibleRows, sizeOfRow);

    populateBoundaryValues(workingArray, sizeOfRow, firstVisibleYIndex,
                           lastVisibleYIndex, numberOfVisibleRows);
    populateBoundaryValues(destinationArray, sizeOfRow, firstVisibleYIndex,
                           lastVisibleYIndex, numberOfVisibleRows);
    if (myRank == 0)
    {
        wholeArray = createArrayOfDoubles(sizeOfRow, sizeOfRow);
        populateBoundaryValues(wholeArray, sizeOfRow, 0,
                               sizeOfRow - 1, sizeOfRow);
    }
    struct RelaxationVariables localRelaxactionVariables;
    localRelaxactionVariables.sizeOfRow = sizeOfRow;
    localRelaxactionVariables.grainSize = grainSize;
    localRelaxactionVariables.precision = userPrecision;
    localRelaxactionVariables.originalArrayPointer = workingArray;
    localRelaxactionVariables.destinationArrayPointer = destinationArray;

    while (countOfPasses == 0)
    {
        performRelaxation(localRelaxactionVariables);
        tempPointer = localRelaxactionVariables.originalArrayPointer;
        localRelaxactionVariables.originalArrayPointer =
            localRelaxactionVariables.destinationArrayPointer;
        localRelaxactionVariables.destinationArrayPointer = tempPointer;
        countOfPasses++;
    }
    //TODO
    /*
    1)  Get one pass then send back to 0 working.
    2)  Get one pass then broadcast if not relaxed, get all other threads
        to stop printing some statement once this is done -- Might not be more efficient. Do it anyway for comparison?
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
    if (myRank == 0)
    {
        int x;
        if (grainSize > 0)
        {
            for (x = 0; x < sizeOfRow * grainSize; x++)
            {
                *(wholeArray[1] + x) =
                    *(localRelaxactionVariables.originalArrayPointer[1] + x);
            }
        }
        /*Noting that MPI_Recv asks for a number of items to receive.
        Importantly it specifies the *maximum* number of items. As such each, if
        carefully programmed so that no other send uses the same tag I can 
        use a larger than actual maximum count here to encapsulate every
        receive from every thread. No data will overlap as each send is careful
        with what data it sends and how much.
        */
        int threadIndex,
            maximumAmount = (grainSize + 1) * sizeOfRow,
            insertYIndex,
            localGrainSize = (sizeOfRow - 2 - remainder) / numberOfThreads;
        ;
        printf("reaminder %d\n", remainder);
        for (threadIndex = 1; threadIndex < nproc; threadIndex++)
        {
            if (remainder)
            {
                //Spreads work across threads, aside from 0 as 0 does more other work.
                if ((remainder - threadIndex) > -1)
                {
                    insertYIndex = threadIndex * (localGrainSize+1);
                }
                else
                {
                    insertYIndex = threadIndex * localGrainSize + 1 + remainder;
                }
            }
            else
            {
                insertYIndex = (threadIndex+1) * grainSize;
            }
            printf("ThreadIndex: %d, InsertIndex: %d\n", threadIndex, insertYIndex);
            MPI_Recv(wholeArray[insertYIndex], maximumAmount,
                     MPI_DOUBLE, threadIndex, threadIndex, MPI_COMM_WORLD,
                     &stat);
        }
    }
    else
    {
        MPI_Send(localRelaxactionVariables.originalArrayPointer[1],
                 grainSize * sizeOfRow, MPI_DOUBLE, 0, myRank, MPI_COMM_WORLD);
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
void performRelaxation(struct RelaxationVariables localVariables)
{
    double **originalArray = localVariables.originalArrayPointer,
           **destinationArray = localVariables.destinationArrayPointer,
           precision = localVariables.precision,
           result;
    int sizeOfRow = localVariables.sizeOfRow,
        grainSize = localVariables.grainSize,
        boundaryXValueIndex = sizeOfRow - 1,
        yIndex, xIndex,
        upperY, lowerY;

    for (yIndex = 1; yIndex <= grainSize; yIndex++)
    {
        upperY = yIndex + 1;
        lowerY = yIndex - 1;
        for (xIndex = 1; xIndex < boundaryXValueIndex; xIndex++)
        {
            result = (originalArray[upperY][xIndex] +
                      originalArray[lowerY][xIndex] +
                      originalArray[yIndex][xIndex + 1] +
                      originalArray[yIndex][xIndex - 1]) /
                     4.0;
            destinationArray[yIndex][xIndex] = result;
        }
    }
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
    int y = 1, x,
        lastVisibleYIndex = grainSize + 1;
    for (y; y < lastVisibleYIndex; y++)
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
void populateBoundaryValues(double **array, int sizeOfRow,
                            int firstVisisbleYIndex,
                            int lastVisibleYIndex, int numberOfVisibleRows)
{
    int y, x, endXIndex = sizeOfRow - 1;

    for (y = 0; y < numberOfVisibleRows; y++)
    {
        array[y][0] = 1.0;
        array[y][endXIndex] = 1.0;
    }
    if (firstVisisbleYIndex == 0)
    {
        for (x = 1; x < endXIndex; x++)
        {
            array[0][x] = 1.0;
        }
    }
    if (lastVisibleYIndex == endXIndex)
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