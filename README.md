# DistributedMemoryParallelism

## Compile

        mpicc -Wall -o parallel ./parallel.c

## Run

        mpirun -np {number of processors} ./parallel

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
* 0
  * Who has sent the data to you.
* 99
  * The tag associated with this message, used to identify messages that are in a stream of other messages.
* MPI_COMM_WORLD
  * What communicator you wish to receive this message through.
* &stat
  * ``` MPI_Status stat ``` A special structure which takes status messages from MPI.