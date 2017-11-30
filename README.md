# DistributedMemoryParallelism

## Compile

        mpicc -Wall -o parallel ./parallel.c

## Run

        mpirun -np {number of processors} ./parallel

-----

### Idea 1

Generate the array in each core, each core works on some segment then
passes the segment back to the main process.

1. Work on segment.
1. Barrier once done.
1. Pass the boundary values to neighbouring processors.
  * Send as one dimensional array?
1. Receive boundary values from neighbouring processors.
1. While this is going on check if is relaxed in processors calculated values.
  * This would need the core to store the previous values in memory.

### Idea 2

Generate the array in the main core, then pass segments to each core. Each
core then passes its segment back.

### Idea 3

One dimensional vs two dimensional?

-----
### MPI Send

``` MPI_Send(n, 1, MPI_INT, 0, 99, MPI_COMM_WORLD); ```

* N
  * Is a pointer to the information you want to send.
* 1
  * The number of items you wish to send.
* MPI_INT
  * What kind of pointer or data you are sending.
* 0
  * Who you are sending it to.
* 99
  * The tag associated with this message, used to identify messages that are in a stream of other messages.
* MPI_COMM_WORLD
  * What communicator you wish to send this through.

### MPI Receive

``` MPI_Recv(n, 1, MPI_INT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &stat); ```

* N
  * Is a pointer to where you want to store the new information.
* 1
  * The number of items you wish to expect to receive.
* MPI_INT
  * What kind of pointer or data you are receiving.
* MPI_ANY_SOURCE
  * Who has sent the data to you, or who you are willing to take it from.
* 99
  * The tag associated with this message, used to identify messages that are in a stream of other messages.
* MPI_COMM_WORLD
  * What communicator you wish to receive this message through.
* &stat
  * ``` MPI_Status stat ``` A special structure which takes status messages from MPI.