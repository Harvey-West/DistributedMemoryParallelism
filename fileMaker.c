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

    for (char threadCount = 1; threadCount <= 48; threadCount += 2)
    {
        for (int arraySize = 1000; arraySize <= 22000; arraySize += 3000)
        {
            for (double precision = 0.1; precision >= 0.001; precision = (precision / (double)10.0))
            {
                char fileName[2048];
                char command[2048];
                sprintf(fileName, "./Files/Thread[%d]Array[%d]Precision[%f].batch", threadCount, arraySize, precision);
                sprintf(command, "mpirun -np %d ../parallel_distributed %d %f", threadCount, arraySize, precision);
                FILE *fp = NULL;
                fp = fopen(fileName, "w");
                if (threadCount > 16)
                {
                    fprintf(fp, "#!/bin/sh\n#SBATCH --account=cm30225\n#SBATCH --partition=teaching\n#SBATCH --nodes=2\n#SBATCH --time=15:00\nmodule load intel/mpi\n");
                }
                else if(threadCount > 32)
                {
                    fprintf(fp, "#!/bin/sh\n#SBATCH --account=cm30225\n#SBATCH --partition=teaching\n#SBATCH --nodes=3\n#SBATCH --time=15:00\nmodule load intel/mpi\n");
                }
                else
                {
                    fprintf(fp, "#!/bin/sh\n#SBATCH --account=cm30225\n#SBATCH --partition=teaching\n#SBATCH --nodes=1\n#SBATCH --time=15:00\nmodule load intel/mpi\n");
                }
                fprintf(fp, command, "a");
            }
        }
    }
    char fileName[2048];
    sprintf(fileName, "./Files/Command.txt");
    FILE *fp = NULL;
    fp = fopen(fileName, "w");
    for (char threadCount = 31; threadCount >= 1; threadCount-=2)
    {
        for (int arraySize = 1000; arraySize <= 22000; arraySize += 3000)
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