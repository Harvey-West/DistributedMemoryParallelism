#include <stdio.h>
#include <mpi.h>
int main(int argc, char **argv)
{
    int rc, myrank, nproc, namelen;
    char name[MPI_MAX_PROCESSOR_NAME];
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        printf("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if (myrank == 0)
    {
        printf("main reports %d procs\n", nproc);
    }
    namelen = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(name, &namelen);
    printf("hello world %d from ’%s’\n", myrank, name);
    /* implicit barrier in Finalize */
    /*MPI_Barrier(MPI_COMM_WORLD);*/
    MPI_Finalize();
    return 0;
}