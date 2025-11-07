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

    // TODO: use MPI_Gather
    long long int* buffer = NULL;
    if (world_rank == 0) {
        buffer = (long long int*)malloc(world_size * sizeof(long long int));
    }
    MPI_Gather(&local_count, 1, MPI_LONG_LONG_INT, buffer, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

    if (world_rank == 0)
    {
        // TODO: PI result
        long long int total_count = 0;
        for (int i = 0; i < world_size; i++) {
            total_count += buffer[i];
        }
        pi_result = 4.0 * (double)total_count / (double)tosses;
        free(buffer);

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
