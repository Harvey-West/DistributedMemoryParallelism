#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <math.h>

struct RelaxParams
{
    int sizeOfArray, startYIndex, noThreads;
    int *isRelaxedPointer;
    double precision;
    double **originalArrayPointer, **destinationArrayPointer;
};
/* A two dimensional array can be expressed by one dimensional malloc whose size
is squared*/
double **createArrayOfDoubles(int sizeOfArray)
{
    int rowIndex;
    double **array = malloc(sizeOfArray * sizeof(double *));
    for (rowIndex = 0; rowIndex < sizeOfArray; rowIndex++)
    {
        array[rowIndex] = malloc(sizeOfArray * sizeof(double));
    }
    return array;
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
/* Fills each index with random double from 0.0 to 1 (not inclusive).
*/
void populateArrayWithRandomNumbers(double **array, int sizeOfArray)
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
/*Performs the operation of taking values from the original array, averages
out 4 adjacent values and sets the same coordinate in the destination array to 
the calculated value.
*/
void relaxArray(struct RelaxParams *args)
{
    double **originalArray = (*args).originalArrayPointer;
    double **destinationArray = (*args).destinationArrayPointer;
    int sizeOfArray = (*args).sizeOfArray;
    int boundaryValueIndex = sizeOfArray - 1;
    double precision = (*args).precision;
    int yIndex, xIndex, startYIndex = (*args).startYIndex,
                        noThreads = (*args).noThreads;
    double result;
    for (yIndex = startYIndex; yIndex < boundaryValueIndex; yIndex += noThreads)
    {
        for (xIndex = 1; xIndex < boundaryValueIndex; xIndex++)
        {
            result = (originalArray[yIndex - 1][xIndex] +
                      originalArray[yIndex + 1][xIndex] +
                      originalArray[yIndex][xIndex + 1] +
                      originalArray[yIndex][xIndex - 1]) /
                     4;
            //If the array is relaxed, calculate if the differenc
            //between the new and old value is less than the precision
            if (*((*args).isRelaxedPointer))
            {
                if (fabs(result - originalArray[yIndex][xIndex]) > precision)
                {
                    *((*args).isRelaxedPointer) = 0;
                }
            }
            destinationArray[yIndex][xIndex] = result;
        }
    }
}
/* For N threads, pass the unique structure of variables to
each thread. Pause the main thread while every pthread computes.
*/
void paralleliseRelaxation(int noThreads,
                           struct RelaxParams **arrayOfRelaxParams,
                           pthread_t *arrayOfThreads)
{
    int threadIndex;

    if (arrayOfThreads == NULL)
    {
        printf("\nERROR: Out of memory creating array of threads.\n");
        exit(1);
    }
    for (threadIndex = 0; threadIndex < noThreads; threadIndex++)
    {
        pthread_create(&arrayOfThreads[threadIndex], NULL,
                       (void *(*)(void *))relaxArray,
                       (void *)&(*arrayOfRelaxParams)[threadIndex]);
    }
    for (threadIndex = 0; threadIndex < noThreads; threadIndex++)
    {
        pthread_join(arrayOfThreads[threadIndex], NULL);
    }
}
/*Creates an array of structures to be passed to each thread. This is used to
dictate what start index each thread will be doing and the needed parameters
for the relaxation to occur.
*/
struct RelaxParams **createThreadRelaxParams(
    struct RelaxParams *originalParams)
{
    int threadIndex, noThreads = (*originalParams).noThreads;
    struct RelaxParams **arrayOfRelaxParams =
        malloc(sizeof(struct RelaxParams *) * noThreads);
    for (threadIndex = 0; threadIndex < noThreads; threadIndex++)
    {
        arrayOfRelaxParams[threadIndex] = malloc(sizeof(struct RelaxParams));
        if (arrayOfRelaxParams[threadIndex] == NULL)
        {
            printf("ERROR: Out of memory generating");
            printf("relaxtion argument pointers");
            exit(1);
        }
        (*arrayOfRelaxParams)[threadIndex].sizeOfArray =
            (*originalParams).sizeOfArray;
        //This is thread index + 1 to exclude already populated boundary
        //values.
        //Further more it ensures each thread has a unique starting index.
        (*arrayOfRelaxParams)[threadIndex].startYIndex = threadIndex + 1;
        (*arrayOfRelaxParams)[threadIndex].isRelaxedPointer =
            (*originalParams).isRelaxedPointer;
        (*arrayOfRelaxParams)[threadIndex].precision =
            (*originalParams).precision;
        (*arrayOfRelaxParams)[threadIndex].originalArrayPointer =
            (*originalParams).originalArrayPointer;
        (*arrayOfRelaxParams)[threadIndex].destinationArrayPointer =
            (*originalParams).destinationArrayPointer;
        (*arrayOfRelaxParams)[threadIndex].noThreads =
            (*originalParams).noThreads;
    }
    return arrayOfRelaxParams;
}
struct RelaxParams *initialiseVariablesForTest(int size, int threads, double prec)
{
    int sizeOfArray = size,
        noThreads = threads,
        isRelaxed = 1;
    double precision = prec;
    double **firstArrayPointer = createArrayOfDoubles(sizeOfArray);
    if (firstArrayPointer == NULL)
    {
        printf("\nERROR: Out of memory creating first array.\n");
        exit(1);
    }

    struct RelaxParams *relaxedParams;
    relaxedParams = malloc(sizeof(struct RelaxParams));
    (*relaxedParams).isRelaxedPointer = malloc(sizeof(int));
    (*relaxedParams).sizeOfArray = sizeOfArray;
    (*relaxedParams).startYIndex = 0;
    (*relaxedParams).noThreads = noThreads;
    *(*relaxedParams).isRelaxedPointer = isRelaxed;
    (*relaxedParams).precision = precision;
    (*relaxedParams).originalArrayPointer = firstArrayPointer;
    return relaxedParams;
}
/*Check input is sane from terminal*/
char isValidInput(int argumentCount, char *userInput[])
{
    char valid = 0;
    if (argumentCount == 4) //isnumber scanf? Put in different function
    {
        //If user provides a string it converts to 0 or 0.0 which is invalid
        //anyway. If user provides 0.1 (for example) followed by a string
        //this is considered valid and will ignore the string after 0.1
        long userSizeOfArray = strtol(userInput[1], NULL, 10),
             userNoThreads = strtol(userInput[2], NULL, 10);
        double userPrecision;
        sscanf(userInput[3], "%lf", &userPrecision);

        // printf("Size of Array: %ld\nNumber of threads: %ld\nPrecision: %0.15f",
        //        userSizeOfArray, userNoThreads, userPrecision);
        // printf("\n\n");
        if (userSizeOfArray < 3 || userSizeOfArray > 20000)
        //Integer limit is much higher, put in more reasonable limit.
        {
            printf("Size of array must be between 3 and 65535.\n");
        }
        else if (userNoThreads < 1 || userNoThreads > 32)
        //Could have more threads but decided on reasonable limit.
        {
            printf("Must have greater than 0 threads and fewer than ");
            printf("33 threads.\n");
        }
        else if (userPrecision >= 1.0 || userPrecision <= 0.00001)
        //Could have a small precision but with larger arrays can take a
        //very long time to run.
        {
            printf("Precision must be between 0.00001 and 1.\n");
        }
        else
        {
            valid = 1;
        }
    }
    else if (argumentCount > 4)
    {
        printf("Too many arguments supplied.\nPlease provide:\nSize of Array");
        printf(":\nNumber of threads\nPrecision.");
    }
    else
    {
        printf("Three arguments expected in order: size of Array, number of t");
        printf("hreads, precision.\n");
    }
    return valid;
}
/* Takes the boundary values from the original array and puts them into the
destination array. This only needs to be done once.
*/
void populateBoundaryIndexes(double **originalArray,
                             double **destinationArray, int sizeOfArray)
{
    int yIndex, xIndex,
        bottomRowIndex = sizeOfArray - 1,
        rightMostColumnIndex = sizeOfArray - 1;
    for (yIndex = 0; yIndex < sizeOfArray; yIndex++)
    {
        //Left column
        destinationArray[yIndex][0] = originalArray[yIndex][0];
        //Right column
        destinationArray[yIndex][rightMostColumnIndex] =
            originalArray[yIndex][rightMostColumnIndex];
    }
    for (xIndex = 0; xIndex < sizeOfArray; xIndex++)
    {
        //Top Row.
        destinationArray[0][xIndex] = originalArray[0][xIndex];
        //Bottom Row.
        destinationArray[bottomRowIndex][xIndex] =
            originalArray[bottomRowIndex][xIndex];
    }
}
int main(int argc, char *argv[])
{
    if (isValidInput(argc, argv) == 1)
    {
        srand(9999);
        int userSizeOfArray = strtol(argv[1], NULL, 10),
            userNoThreads = strtol(argv[2], NULL, 10),
            count = 1;
        double userPrecision, timeSpent;
        double **tempPointer;
        //Structure to calculate relaxation time by storing time
        //spent in main function.
        struct timespec startTime, endTime;
        sscanf(argv[3], "%lf", &userPrecision);

        struct RelaxParams *relaxedParams = initialiseVariablesForTest(
            userSizeOfArray, userNoThreads, userPrecision);

        populateArrayWithRandomNumbers(
            (*relaxedParams).originalArrayPointer,
            (*relaxedParams).sizeOfArray);

        printArray((*relaxedParams).originalArrayPointer,
                   (*relaxedParams).sizeOfArray);

        //Start timing my algoritm.
        clock_gettime(CLOCK_MONOTONIC, &startTime);
        //Create second array, the destination array.
        (*relaxedParams).destinationArrayPointer =
            createArrayOfDoubles(userSizeOfArray);
        if ((*relaxedParams).destinationArrayPointer == NULL)
        {
            printf("\nERROR: Out of memory creating second array.\n");
            exit(1);
        }
        //Create an array of unique structures to pass to each array.
        //Only needs to be called once.
        struct RelaxParams **arrayOfRelaxParams =
            createThreadRelaxParams(relaxedParams);
        
        //Create an array of pthread pointers.
        pthread_t *arrayOfThreads = malloc(sizeof(pthread_t) * userNoThreads);

        //Only needs to be called once.
        populateBoundaryIndexes((*relaxedParams).originalArrayPointer,
                                (*relaxedParams).destinationArrayPointer,
                                (*relaxedParams).sizeOfArray);
        /************/
        //Perform the first pass.
        paralleliseRelaxation((*relaxedParams).noThreads,
                              arrayOfRelaxParams,
                              arrayOfThreads);

        while (*(*relaxedParams).isRelaxedPointer == 0)
        {
            //Reset isRelaxed to check array after next relaxation.
            *(*relaxedParams).isRelaxedPointer = 1;

            //Switch original and destination array pointers.
            //This allows for only 2 mallocs, one for each array.
            tempPointer = (*relaxedParams).originalArrayPointer;
            (*relaxedParams).originalArrayPointer =
                (*relaxedParams).destinationArrayPointer;
            (*relaxedParams).destinationArrayPointer = tempPointer;

            int threadIndex;
            //Swap pointers around for each structure
            for (threadIndex = 0;
                 threadIndex < (*relaxedParams).noThreads; threadIndex++)
            {
                (*arrayOfRelaxParams)[threadIndex].destinationArrayPointer =
                    (*relaxedParams).destinationArrayPointer;
                (*arrayOfRelaxParams)[threadIndex].originalArrayPointer =
                    (*relaxedParams).originalArrayPointer;
            }

            //Perform relaxation with arrays swapped.
            paralleliseRelaxation((*relaxedParams).noThreads,
                                  arrayOfRelaxParams,
                                  arrayOfThreads);
            count++;
        }

        clock_gettime(CLOCK_MONOTONIC, &endTime);
        timeSpent = endTime.tv_sec - startTime.tv_sec;
        timeSpent += (endTime.tv_nsec - startTime.tv_nsec) / 1000000000.0;
        if (count % 2)
        {
            printArray((*relaxedParams).originalArrayPointer,
                       (*relaxedParams).sizeOfArray);
        }
        else
        {
            printArray((*relaxedParams).destinationArrayPointer,
                       (*relaxedParams).sizeOfArray);
        }
        // printf("%d, %d, %d, %f\n", userSizeOfArray, userNoThreads, count, timeSpent);
        printf("Parallel\nRow Size:%d, Iteration Count:%d, Time Spent:%f, Thread Count:%d\n", userSizeOfArray, count, timeSpent, userNoThreads);
        // printf("%f,",timeSpent);
    }
    return 0;
}