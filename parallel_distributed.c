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
        on, global index 5.

        Thread 1 must be able to see the last row that thread 0 works on, i.e.
        the first visible Y from the full array that thread 1 sees. This
        corresponds to thread index in the full array of 4.
        Thread 1 will need to see 4, to work on row indexes 5.
        Thread 1 will also need to see the final boundary row 9, however will
        not work on this row as it never changes like row 0, to calculate
        global row index 8. This is not included in the array as it is just:

            First visible Y + local grain size + 1

        The resultant array for this example would be:

        Array[0][0] - Grain Size for Thread 0 = 4
        Array[0][1] - First Visible Y index for Thread 0 = 0

        Array[1][0] - Grain Size for Thread 1 = 4
        Array[1][1] - First Visible Y index for Thread 1 = 4

        Justification for this is it makes it easier to:
            - each thread to swap boundary values
            - final reduction step to know where to insert
                rows into the final result array.
            - when generating local array, know if there are boundary values
                to generate.
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
        //row past the final grain index.
        numberOfVisibleRows = localGrainSize + 2;

    double **workingArray, **destinationArray, //These arrays are for performing
    //the relaxation.
     **wholeArray,//End result will be put into the "wholeArray" on whichever
     //core calls createArrayOfDoubles on this variable. 
     **tempPointer//Used to swap pointers in the structure passed to function
     //which relaxes values from one array and puts result into the other. Means
     //that malloc only needs to be called once.
     ;

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
        //If not pseudo main, send "if relaxed" to pseudo main.
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
        //If any of the threads say it is not relaxed, broadcast to all that
        //process of relaxing has to be repeated.
        if (numberOfThreads > 1)
        {
            MPI_Bcast(&globalIsRelaxed, 1, MPI_INT,
                      lastThreadRank, MPI_COMM_WORLD);
        }

        //Ensure every thread has received the broadcast
        MPI_Barrier(MPI_COMM_WORLD);
        /*
        If not relaxed, and more than one thread, then send and receive edge
        rows to nearest neighbours.
        
        A sending edge row is the first and last grain that the local thread
        worked on.

        A receiving edge row is not within the grain, and was not relaxed
        in the local thread, but the boundaries around this grain.
        */
        if (!globalIsRelaxed && numberOfThreads > 1)
        {
            /*
                Non blocking send and blocking receive ensures that there is
                no halting problem where two cores are waiting for the other to
                send before sending themselves. Only works if the receive is
                blocking because a blocking send and non blocking receive 
                could mean that the data is not received before another
                iteration of performRelaxation occurs. 
            */
            int tagOffset = (numberOfThreads + 1);
            if (myRank == 0)
            {
                //If thread 0, then there is only one overlapping row. Send the
                //last grain and receive for the last grain index plus one.
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
                //If local thread has the greatest thread rank, then send
                //the first working grain index and receive in the first index
                //at 0.
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
                //See above commends, sends edge rows to nearest neighbour.
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
                         tagOffset + firstVisibleYIndex +
                             numberOfVisibleRows - 1,
                         MPI_COMM_WORLD, &stat);
                MPI_Recv(localRelaxactionVariables.originalArrayPointer[0],
                         sizeOfRow, MPI_DOUBLE, myRank - 1,
                         tagOffset + firstVisibleYIndex, MPI_COMM_WORLD, &stat);
            }
        }
    }

    //While loop has exited, as such relaxation is complete.
    //Results need to be collected as final stage.

    //If local thread is the pseudo main, receive all of the results from
    //other threads and collate them.
    if (numberOfThreads > 1)
    {
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
    }
    else
    {
        int x;
        for (x = 0; x < sizeOfRow * localGrainSize; x++)
        {
            *(wholeArray[firstVisibleYIndex + 1] + x) =
                *(localRelaxactionVariables.originalArrayPointer[1] + x);
        }
    }

    //Ensure all messages are sent and that the pseudo main thread has received
    //all results
    if (myRank == lastThreadRank)
    {
        clock_gettime(CLOCK_MONOTONIC, &endTime);
        timeSpent = endTime.tv_sec - startTime.tv_sec;
        timeSpent += (endTime.tv_nsec - startTime.tv_nsec) / 1000000000.0;
        //printf("\nNumber of passes: %d, TimeSpent: %f\n", countOfPasses, timeSpent);
        printf("%d, %d, %d, %f\n\n", sizeOfRow, numberOfThreads, countOfPasses, timeSpent);
        // printArray(wholeArray, sizeOfRow);
        // printGrain(workingArray, myRank, grainSize, sizeOfRow, numberOfThreads);
        // printVisibleArea(workingArray, numberOfVisibleRows, sizeOfRow);
    }
    if (numberOfThreads > 1)
    {
        MPI_Finalize();
    }
    return 0;
}
/*
Given the structure of type RelaxationVariables, take values from one array,
relax and place them in the destination array. Check if the new value is relaxed
or not.
Once relaxation iteration is complete, return if grain is relaxed compared to
the precision or not.
*/
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
        //Minor optimisation means these two values are calculated fewer times
        //than if they were included in the following for loop.
        upperY = yIndex + 1;
        lowerY = yIndex - 1;
        for (xIndex = 1; xIndex < boundaryXValueIndex; xIndex++)
        {
            //Get relevant values and divide by 4.0 to relax.
            result = (originalArray[upperY][xIndex] +
                      originalArray[lowerY][xIndex] +
                      originalArray[yIndex][xIndex + 1] +
                      originalArray[yIndex][xIndex - 1]) /
                     4.0;
            //Optimisation that means calculating if relaxed may not need to
            //occur. In the worst case every value needs to both check
            //isRelaxed and then calculate if it is relaxed. However the cost
            //is less than always calculating and never checking if already
            //found a value that is relaxed.
            //Furthermore spreading this computation amongst cores improved
            //speedup, rather than collecting and having one core calculate
            //if the entire array was relaxed. More CPUs spent time idle and
            //the program was slower, this method alleviates that issue.
            if (isRelaxed)
            {
                //Check if result changes more than the precision.
                if (fabs(result - originalArray[yIndex][xIndex]) > precision)
                {
                    isRelaxed = 0;
                }
            }
            //Populate destination array with results.
            destinationArray[yIndex][xIndex] = result;
        }
    }
    //Return if local grain is relaxed or not. The result values are already
    //assigned via use of the pointers.
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
Two dimensional access is useful because it reduces the required computations
to determine if the value at an index is a boundary or not.
*/
double **createArrayOfDoubles(int numberOfRows, int sizeOfRow)
{
    int rowIndex;
    //Initialise all values to 0.0
    double *array = calloc(numberOfRows * sizeOfRow, sizeof(double *));
    double **rowPointers = malloc(sizeOfRow * sizeof(double *));
    for (rowIndex = 0; rowIndex < numberOfRows; rowIndex++)
    {
        rowPointers[rowIndex] = &array[rowIndex * sizeOfRow];
    }
    return rowPointers;
}

/* Fills each boundary index with 1.0. A fixed value is used to ensure
accurate comparison between array sizes.*/
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
    //If the grain can not be split equally amongst threads, segment the work
    //between all threads. The maximum remainder can only be numberOfThreads-1.
    //Due to this effect, the thread with the greatest rank in MPI_COMM_WORLD
    //will be the "pseudo main" in that they will perform any additional work to
    //compensate for potentially less work.
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
//Receives results from all other threads and collates them into wholeArray.
void receiveFromAll(int sizeOfRow, int numberOfThreads, double ***wholeArray,
                    int (*threadGrainAndFirstVisibleYIndexArray)[2])
{
    MPI_Status stat;
    int threadIndex,
        insertYIndex;
    for (threadIndex = 0; threadIndex < numberOfThreads - 1; threadIndex++)
    {
        insertYIndex = threadGrainAndFirstVisibleYIndexArray
                           [threadIndex][1] +
                       1;
        MPI_Recv((*wholeArray)[insertYIndex],
                 threadGrainAndFirstVisibleYIndexArray
                         [threadIndex][0] *
                     sizeOfRow,
                 MPI_DOUBLE, threadIndex, threadIndex, MPI_COMM_WORLD,
                 &stat);
    }
}

/*
This function returns what Y index, from a global view of the whole array, that
a thread will see. Each thread will be working on some number of rows:

1.0, 1.0, 1.0
0.0, 0.0, 0.0
1.0, 1.0, 1.0

In this case if the grain size is 1, i.e. it is only going to relax the second
row of 0.0. It will need to see the row preceding it, the first 1.0 row.
This function returns what Y index the first 1.0 row would have if it was taken
from the entire result array.
*/
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