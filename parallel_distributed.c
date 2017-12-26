#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
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
                    int (*threadGrainAndFirstVisibleYIndexArray)[2]);
int getFirstVisibleYIndex(int remainder, int myRank, int grainSize);

int main(int argc, char **argv)
{
    //Start - Initialise variables for test
    int sizeOfRow = strtol(argv[1], NULL, 10),
        //Each thread has an amount of work it will do, this is the grain.
        //For my program the grain size is the number of rows that the thread
        //will relax.
        localGrainSize, remainder, threadIndex, countOfPasses = 0,
        rc, myRank, numberOfThreads, lastThreadRank,
        globalIsRelaxed = 0, localIsRelaxed = 1;
    double userPrecision, timeSpent;
    struct timespec startTime, endTime;
    sscanf(argv[2], "%lf", &userPrecision);

    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        printf("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Status stat;
    MPI_Request req;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfThreads);

    remainder = (sizeOfRow - 2) % numberOfThreads;

    /*The last thread rank will be used a psuedo main.
    What I mean by this is that, the last thread will handle additional
    computation beyond relaxing the given grain.
    If there is a remainder, this lastThreadRank will have a smaller grain size
    compared to its peers to try and accomodate this "extra" pseudo main work.
    */
    lastThreadRank = numberOfThreads - 1;

    /*For each thread, create an array that will contain:
    Each index through the array corresponds to the thread rank.
        - Grain Size: The number of rows to relax
        - First Visible Y Index: For the sake of knowing where each thread 
        starts and what rows they need to see from a global view.
        
        If we say that the array size is 10x10 with 2 threads trying
        to relax it. We know that the first and last row is the boundary and
        will never change. As such the possible working area is actually 8 rows
        and 10 columns (NB could do 8 columns and 10 rows but not done in this
        case, same principle).

        Both threads get an equal grain size of 4 (8 / 2).

        Thread 0 must be able to see the original y index 0 so that they can
        relax rows 1,2,3,4.
        Additionally, thread 0 needs to see their last grain row index + 1 to 
        calculate row 4. This corresponds to the first row that thread 1 works
        on, 5.

        Thread 1 must be able to see the last row that thread 0 works on, i.e.
        the first visible Y from the full array that thread 1 sees is 4.
        Thread 1 will need to see 4, to work on row indexes 5,6,7,8.
        Thread 1 will also need to see the final boundary row 9, however will
        not work on this row as it never changes.

        The resultant array for this example would be:

        Array[0][0] - Grain Size for Thread 0 = 4
        Array[0][1] - First Visible Y index for Thread 0 = 0

        Array[1][0] - Grain Size for Thread 1 = 4
        Array[1][1] - First Visible Y index for Thread 1 = 4
    */
    int threadGrainAndFirstVisibleYIndexArray[numberOfThreads][2];

    for (threadIndex = 0; threadIndex < numberOfThreads; threadIndex++)
    {
        threadGrainAndFirstVisibleYIndexArray[threadIndex][0] =
            getGrainSize(sizeOfRow, numberOfThreads, threadIndex);
        threadGrainAndFirstVisibleYIndexArray[threadIndex][1] =
            getFirstVisibleYIndex(
                remainder, threadIndex,
                threadGrainAndFirstVisibleYIndexArray[threadIndex][0]);
    }

    //We have a list of every other thread's grain and first visible,
    //get the local to this thread values for these.
    localGrainSize = threadGrainAndFirstVisibleYIndexArray[myRank][0];
    int firstVisibleYIndex = threadGrainAndFirstVisibleYIndexArray[myRank][1],
        //As outlined above, the number of total visible rows will be the grain
        //size plus one for the row above and plus one additional one for the
        //row belowthe grain.
        numberOfVisibleRows = localGrainSize + 2;

    double **workingArray, **destinationArray, **wholeArray, **tempPointer;

    //These arrays are malloced once with values taken from one and results
    //put into the other.
    workingArray = createArrayOfDoubles(numberOfVisibleRows, sizeOfRow);
    destinationArray = createArrayOfDoubles(numberOfVisibleRows, sizeOfRow);

    //Fill values of array. Using the first visible Y we can determine if they
    //are boundary values or not and as such fill them in.
    populateBoundaryValues(workingArray, sizeOfRow, firstVisibleYIndex,
                           numberOfVisibleRows);
    populateBoundaryValues(destinationArray, sizeOfRow, firstVisibleYIndex,
                           numberOfVisibleRows);

    //The psuedo main thread will collect all the results into wholeArray and
    //time the test.
    if (myRank == lastThreadRank)
    {
        wholeArray = createArrayOfDoubles(sizeOfRow, sizeOfRow);
        populateBoundaryValues(wholeArray, sizeOfRow, 0, sizeOfRow);
        clock_gettime(CLOCK_MONOTONIC, &startTime);
    }

    struct RelaxationVariables localRelaxactionVariables;
    localRelaxactionVariables.sizeOfRow = sizeOfRow;
    localRelaxactionVariables.grainSize = localGrainSize;
    localRelaxactionVariables.precision = userPrecision;
    localRelaxactionVariables.originalArrayPointer = workingArray;
    localRelaxactionVariables.destinationArrayPointer = destinationArray;

    //While all threads say that the array is not relaxed.
    while (!globalIsRelaxed)
    {
        localIsRelaxed = performRelaxation(localRelaxactionVariables);

        //Swap the pointers around for the next relaxation.
        tempPointer = localRelaxactionVariables.originalArrayPointer;
        localRelaxactionVariables.originalArrayPointer =
            localRelaxactionVariables.destinationArrayPointer;
        localRelaxactionVariables.destinationArrayPointer = tempPointer;
        countOfPasses++;
        if (myRank < lastThreadRank)
        {
            MPI_Send(&localIsRelaxed, 1, MPI_INT, lastThreadRank,
                     numberOfThreads, MPI_COMM_WORLD);
        }
        else
        {
            globalIsRelaxed = localIsRelaxed;
            for (threadIndex = 0; threadIndex < lastThreadRank; threadIndex++)
            {
                MPI_Recv(&localIsRelaxed, 1, MPI_INT, threadIndex,
                         numberOfThreads, MPI_COMM_WORLD, &stat);
                if (localIsRelaxed == 0)
                {
                    globalIsRelaxed = 0;
                }
            }
        }
        MPI_Bcast(&globalIsRelaxed, 1, MPI_INT, lastThreadRank, MPI_COMM_WORLD);
        if (!globalIsRelaxed && numberOfThreads > 1)
        {
            int tagOffset = (numberOfThreads + 1);
            if (myRank == 0)
            {
                MPI_Isend(localRelaxactionVariables.originalArrayPointer
                              [localGrainSize],
                          sizeOfRow, MPI_DOUBLE, myRank + 1,
                          tagOffset + firstVisibleYIndex + localGrainSize,
                          MPI_COMM_WORLD, &req);
                MPI_Recv(localRelaxactionVariables.originalArrayPointer
                             [numberOfVisibleRows - 1],
                         sizeOfRow, MPI_DOUBLE, myRank + 1,
                         (tagOffset + firstVisibleYIndex +
                          numberOfVisibleRows - 1),
                         MPI_COMM_WORLD, &stat);
            }
            else if (myRank == lastThreadRank)
            {
                MPI_Isend(localRelaxactionVariables.originalArrayPointer[1],
                          sizeOfRow, MPI_DOUBLE, myRank - 1,
                          tagOffset + firstVisibleYIndex + 1,
                          MPI_COMM_WORLD, &req);
                MPI_Recv(localRelaxactionVariables.originalArrayPointer[0],
                         sizeOfRow, MPI_DOUBLE, myRank - 1,
                         tagOffset + firstVisibleYIndex, MPI_COMM_WORLD, &stat);
            }
            else
            {

                MPI_Isend(localRelaxactionVariables.originalArrayPointer[1],
                          sizeOfRow, MPI_DOUBLE, myRank - 1,
                          tagOffset + firstVisibleYIndex + 1,
                          MPI_COMM_WORLD, &req);
                MPI_Isend(localRelaxactionVariables.originalArrayPointer
                              [localGrainSize],
                          sizeOfRow, MPI_DOUBLE, myRank + 1,
                          tagOffset + firstVisibleYIndex + localGrainSize,
                          MPI_COMM_WORLD, &req);
                MPI_Recv(localRelaxactionVariables.originalArrayPointer
                             [numberOfVisibleRows - 1],
                         sizeOfRow, MPI_DOUBLE, myRank + 1,
                         tagOffset + firstVisibleYIndex + numberOfVisibleRows - 1, MPI_COMM_WORLD, &stat);
                MPI_Recv(localRelaxactionVariables.originalArrayPointer[0],
                         sizeOfRow, MPI_DOUBLE, myRank - 1,
                         tagOffset + firstVisibleYIndex, MPI_COMM_WORLD, &stat);
            }
        }
    }

    if (myRank == lastThreadRank)
    {

        int x;
        for (x = 0; x < sizeOfRow * localGrainSize; x++)
        {
            *(wholeArray[firstVisibleYIndex + 1] + x) =
                *(localRelaxactionVariables.originalArrayPointer[1] + x);
        }
        receiveFromAll(sizeOfRow, numberOfThreads, &wholeArray,
                       (int(*)[2]) & threadGrainAndFirstVisibleYIndexArray);
    }
    else
    {
        MPI_Send(localRelaxactionVariables.originalArrayPointer[1],
                 localGrainSize * sizeOfRow, MPI_DOUBLE, lastThreadRank,
                 myRank, MPI_COMM_WORLD);
    }
    /* implicit barrier in Finalize */
    /*MPI_Barrier(MPI_COMM_WORLD);*/
    MPI_Finalize();
    if (myRank == lastThreadRank)
    {
        clock_gettime(CLOCK_MONOTONIC, &endTime);
        timeSpent = endTime.tv_sec - startTime.tv_sec;
        timeSpent += (endTime.tv_nsec - startTime.tv_nsec) / 1000000000.0;
        printf("\nNumber of passes: %d, TimeSpent: %f\n", countOfPasses, timeSpent);
        // printArray(wholeArray, sizeOfRow);
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
    for (y = 1; y < lastVisibleYIndex; y++)
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
        if ((remaining - myRank) > 0)
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
                    int (*threadGrainAndFirstVisibleYIndexArray)[2])
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
        maximumAmount = (threadGrainAndFirstVisibleYIndexArray
                             [numberOfThreads - 1][0] +
                         1) *
                        sizeOfRow,
        insertYIndex;
    for (threadIndex = 0; threadIndex < numberOfThreads - 1; threadIndex++)
    {
        insertYIndex = threadGrainAndFirstVisibleYIndexArray[threadIndex][1] + 1;
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
        if (myRank >= remainder)
        {
            firstVisibleYIndex = myRank * grainSize + remainder;
        }
        else
        {
            firstVisibleYIndex = myRank * grainSize;
        }
    }
    else
    {
        firstVisibleYIndex = myRank * grainSize;
    }
    return firstVisibleYIndex;
}