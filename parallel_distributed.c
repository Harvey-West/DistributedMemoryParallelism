#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
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
struct RelaxationVariables
{
    int sizeOfRow, grainSize;
    double precision;
    double **originalArrayPointer, **destinationArrayPointer;
};

double **createArrayOfDoubles(int numberOfRows, int sizeOfRow);
void populateBoundaryValues(double **array, int sizeOfRow, int startYIndex,
                            int numberOfVisibleRows);
int getGrainSize(int sizeOfRow, int numberOfThreads, int myRank);
int performRelaxation(struct RelaxationVariables);
void printArray(double **array, int sizeOfRow);
void printGrain(double **array, int offset, int grainSize, int sizeOfRow,
                int numberOfThreads);
void printVisibleArea(double **array, int numberOfVisibleRows, int sizeOfRow);
void receiveFromAll(int sizeOfRow, int numberOfThreads, double ***wholeArray,
                    int threadGrainAndFirstVisibleYIndexArray[][2]);
int getFirstVisibleYIndex(int remainder, int myRank, int grainSize);
int main(int argc, char **argv)
{
    int sizeOfRow = strtol(argv[1], NULL, 10),
        localGrainSize,
        remainder,
        countOfPasses = 0,
        rc, myRank, numberOfThreads,
        globalIsRelaxed = 0, localIsRelaxed = 1;
    double userPrecision;
    sscanf(argv[2], "%lf", &userPrecision);

    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        printf("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Status stat;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfThreads);
    remainder = (sizeOfRow - 2) % numberOfThreads;
    int threadIndex, threadGrainAndFirstVisibleYIndexArray[numberOfThreads][2];
    for (threadIndex = 0; threadIndex < numberOfThreads; threadIndex++)
    {
        threadGrainAndFirstVisibleYIndexArray[threadIndex][0] =
            getGrainSize(sizeOfRow, numberOfThreads, threadIndex);
        threadGrainAndFirstVisibleYIndexArray[threadIndex][1] =
            getFirstVisibleYIndex(
                remainder, threadIndex,
                threadGrainAndFirstVisibleYIndexArray[threadIndex][0]);
    }
    localGrainSize = threadGrainAndFirstVisibleYIndexArray[myRank][0];
    int firstVisibleYIndex = threadGrainAndFirstVisibleYIndexArray[myRank][1],
        numberOfVisibleRows = localGrainSize + 2;

    double **workingArray, **destinationArray, **wholeArray, **tempPointer;

    workingArray = createArrayOfDoubles(numberOfVisibleRows, sizeOfRow);
    destinationArray = createArrayOfDoubles(numberOfVisibleRows, sizeOfRow);

    populateBoundaryValues(workingArray, sizeOfRow, firstVisibleYIndex,
                           numberOfVisibleRows);
    populateBoundaryValues(destinationArray, sizeOfRow, firstVisibleYIndex,
                           numberOfVisibleRows);

    if (myRank == 0)
    {
        wholeArray = createArrayOfDoubles(sizeOfRow, sizeOfRow);
        populateBoundaryValues(wholeArray, sizeOfRow, 0, sizeOfRow);
    }

    struct RelaxationVariables localRelaxactionVariables;
    localRelaxactionVariables.sizeOfRow = sizeOfRow;
    localRelaxactionVariables.grainSize = localGrainSize;
    localRelaxactionVariables.precision = userPrecision;
    localRelaxactionVariables.originalArrayPointer = workingArray;
    localRelaxactionVariables.destinationArrayPointer = destinationArray;

    //while(globalIsRelaxed)
    while (!globalIsRelaxed)
    {
        localIsRelaxed = performRelaxation(localRelaxactionVariables);
        tempPointer = localRelaxactionVariables.originalArrayPointer;
        localRelaxactionVariables.originalArrayPointer =
            localRelaxactionVariables.destinationArrayPointer;
        localRelaxactionVariables.destinationArrayPointer = tempPointer;
        countOfPasses++;
        if (myRank)
        {
            MPI_Send(&localIsRelaxed, 1, MPI_INT, 0, 555, MPI_COMM_WORLD);
        }
        else
        {
            globalIsRelaxed = localIsRelaxed;
            for (threadIndex = 1; threadIndex < numberOfThreads; threadIndex++)
            {
                MPI_Recv(&localIsRelaxed, 1, MPI_INT, threadIndex, 555,
                         MPI_COMM_WORLD, &stat);
                if (localIsRelaxed == 1)
                {
                    globalIsRelaxed = 1;
                }
            }
        }
        MPI_Bcast(&globalIsRelaxed, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(!globalIsRelaxed)
        {
            //Send to nearest neighbour grain index 0 and grainsize +1;
        }
    }

    printf("hello world from ’%d’, grainSize = %d, startY = %d\n", myRank, localGrainSize, firstVisibleYIndex);
    MPI_Barrier(MPI_COMM_WORLD);

    if (myRank == 0)
    {
        //In some cases where the offset of trying to prevent thread 0 from
        //getting too much work, it may have no grain to work on.
        if (localGrainSize > 0)
        {
            int x;
            for (x = 0; x < sizeOfRow * localGrainSize; x++)
            {
                *(wholeArray[1] + x) =
                    *(localRelaxactionVariables.originalArrayPointer[1] + x);
            }
        }
        receiveFromAll(sizeOfRow, numberOfThreads, &wholeArray,
                       &threadGrainAndFirstVisibleYIndexArray);
    }
    else
    {
        MPI_Send(localRelaxactionVariables.originalArrayPointer[1],
                 localGrainSize * sizeOfRow, MPI_DOUBLE, 0,
                 myRank, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* implicit barrier in Finalize */
    /*MPI_Barrier(MPI_COMM_WORLD);*/
    MPI_Finalize();
    if (myRank == 0)
    {
        printArray(wholeArray, sizeOfRow);
        // printGrain(workingArray, myRank, grainSize, sizeOfRow, numberOfThreads);
        // printVisibleArea(workingArray, numberOfVisibleRows, sizeOfRow);
    }
    return 0;
}
int performRelaxation(struct RelaxationVariables localVariables)
{
    double **originalArray = localVariables.originalArrayPointer,
           **destinationArray = localVariables.destinationArrayPointer,
           precision = localVariables.precision,
           result;
    int sizeOfRow = localVariables.sizeOfRow,
        grainSize = localVariables.grainSize,
        boundaryXValueIndex = sizeOfRow - 1,
        yIndex, xIndex,
        upperY, lowerY,
        isRelaxed = 1;
    MPI_Status stat;
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
            if (isRelaxed)
            {
                if (fabs(result - originalArray[yIndex][xIndex]) > precision)
                {
                    isRelaxed = 0;
                }
            }
            destinationArray[yIndex][xIndex] = result;
        }
    }
    return isRelaxed;
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
                int numberOfThreads)
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

/* A two dimensional array can be expressed by one dimensional malloc whose size
is squared.
To enable two dimensional access, simply malloc pointers to each row index.
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

/* Fills each boundary index with 1.0*/
void populateBoundaryValues(double **array, int sizeOfRow,
                            int firstVisisbleYIndex, int numberOfVisibleRows)
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
    if (firstVisisbleYIndex + numberOfVisibleRows == sizeOfRow)
    {
        for (x = 1; x < endXIndex; x++)
        {
            array[numberOfVisibleRows - 1][x] = 1.0;
        }
    }
}

/*Returns the number of rows that a thread will produce results for,
i.e. the grain.*/
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
void receiveFromAll(int sizeOfRow, int numberOfThreads, double ***wholeArray,
                    int threadGrainAndFirstVisibleYIndexArray[][2])
{
    MPI_Status stat;
    /*Noting that MPI_Recv asks for a number of items to receive.
        Importantly it specifies the *maximum* number of items. As such each, if
        carefully programmed so that no other send uses the same tag I can 
        use a larger than actual maximum count here to encapsulate every
        receive from every thread. No data will overlap as each send is careful
        with what data it sends and how much.
        */
    int threadIndex,
        maximumAmount = (threadGrainAndFirstVisibleYIndexArray[1][0] + 1) *
                        sizeOfRow,
        insertYIndex;
    for (threadIndex = 1; threadIndex < numberOfThreads; threadIndex++)
    {
        insertYIndex = threadGrainAndFirstVisibleYIndexArray[threadIndex][1] + 1;
        printf("ThreadIndex: %d, InsertIndex: %d\n", threadIndex, insertYIndex);
        MPI_Recv((*wholeArray)[insertYIndex], maximumAmount,
                 MPI_DOUBLE, threadIndex, threadIndex, MPI_COMM_WORLD,
                 &stat);
    }
}

int getFirstVisibleYIndex(int remainder, int myRank, int grainSize)
{
    int firstVisibleYIndex = 0;
    if (remainder)
    {
        if (myRank > remainder)
        {
            firstVisibleYIndex = myRank * grainSize + remainder;
        }
        else if (myRank == 0)
        {
            firstVisibleYIndex = myRank * grainSize;
        }
        else
        {
            firstVisibleYIndex = myRank * grainSize - 1;
        }
    }
    else
    {
        firstVisibleYIndex = myRank * grainSize;
    }
    return firstVisibleYIndex;
}