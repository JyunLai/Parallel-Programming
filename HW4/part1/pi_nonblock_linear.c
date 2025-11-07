#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---

    // TODO: MPI init
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    long long int my_tosses = tosses / world_size;
    long long int remainder = tosses / world_size;
    if (world_rank < remainder) {
        my_tosses++;
    }
    
    unsigned int seed = time(NULL) + world_rank;

    long long int local_count = 0;
    for (long long int i = 0; i < my_tosses; i++) {
        double x = (double)rand_r(&seed) / RAND_MAX;
        double y = (double)rand_r(&seed) / RAND_MAX;
        if ((x*x + y*y) <= 1.0) {
            local_count++;
        }
    }
    long long int total_count = 0;

    if (world_rank > 0)
    {
        // TODO: MPI workers
        MPI_Send(&local_count, 1, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
    }
    else if (world_rank == 0)
    {
        // TODO: non-blocking MPI communication.
        total_count = local_count;
        int num_worker = world_size - 1;

        if (num_worker > 0) {
            MPI_Request* request = (MPI_Request*)malloc(num_worker * sizeof(MPI_Request));
            MPI_Status* status = (MPI_Status*)malloc(num_worker * sizeof(MPI_Status));
            long long int* received_count = (long long int*)malloc(num_worker * sizeof(long long int));

            for (int i = 0; i < num_worker; i++) {
                int source_rank = i + 1;
                MPI_Irecv(&received_count[i], 1, MPI_LONG_LONG_INT, source_rank, 0, MPI_COMM_WORLD, &request[i]);
            }
            MPI_Waitall(num_worker, request, status);
            for (int i = 0; i < num_worker; i++) {
                total_count += received_count[i];
            }
            free(request);
            free(status);
            free(received_count);
        }
    }

    if (world_rank == 0)
    {
        // TODO: PI result
        pi_result = 4.0 * (double)total_count / (double)tosses;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
