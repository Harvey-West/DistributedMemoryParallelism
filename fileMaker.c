#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #!/bin/sh
// # Account & partition (must have these)
// #SBATCH --account=cm30225
// #SBATCH --partition=teaching
// # Name of job (optional)
// #SBATCH --job-name=2K0
// # one node
// #SBATCH --nodes=1
// # Run the program
int main(int argc, char *argv[])
{

    for (char threadCount = 1; threadCount <= 16; threadCount++)
    {
        for (int arraySize = 1000; arraySize <= 21000; arraySize += 2000)
        {
            for (double precision = 0.1; precision >= 0.001; precision = (precision / (double) 10.0))
            {
                char fileName[2048];
                char command[2048];
                sprintf(fileName, "./Files/Thread[%d]Array[%d]Precision[%f].batch", threadCount, arraySize, precision);
                sprintf(command, "mpirun -np %d ../parallel_distributed %d %f", threadCount, arraySize, precision);
                FILE *fp = NULL;
                fp = fopen(fileName, "w");
                fprintf(fp, "#!/bin/sh\n#SBATCH --account=cm30225\n#SBATCH --partition=teaching\n#SBATCH --nodes=1\n#SBATCH --time=15:00\n");
                fprintf(fp, command, "a");
            }
        }
    }
    char fileName[2048];
    sprintf(fileName, "./Files/Command.txt");
    FILE *fp = NULL;
    fp = fopen(fileName, "w");
    for (char threadCount = 16; threadCount >= 1; threadCount--)
    {
        for (int arraySize = 1000; arraySize <= 21000; arraySize += 2000)
        {
            for (double precision = 0.1; precision >= 0.001; precision = (precision / 10.0))
            {
                char command[2048];
                sprintf(command, "sbatch ./Thread[%d]Array[%d]Precision[%f].batch\n", threadCount, arraySize, precision);
                fprintf(fp, command, "a");
            }
        }
    }
}