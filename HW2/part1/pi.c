#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "montecarlo.h"

int main(int argc, char **argv){
    if (argc != 3){
        fprintf(stderr, "usage: %s <num_threads> <num_tosses>\n", argv[0]);
        return 1;
    }

    int threads = atoi(argv[1]);
    long long tosses = atoll(argv[2]);
    if (threads <= 0 || tosses <= 0) return 1;

    omp_set_num_threads(threads);

    long long in_circle = count_in_circle_parallel(tosses);

    double pi = 4.0 * (double)in_circle / (double)tosses;
    printf("%lf\n", pi);

    return 0;
}