#include "montecarlo.h"
#include <stdlib.h>
#include <omp.h>
#include <time.h>

long long count_in_circle_parallel(long long tosses){
    long long in_circle = 0;
    #pragma omp parallel
    {
        unsigned int lcg_state = (unsigned int)time(NULL) ^ ((unsigned int)omp_get_thread_num() << 16);
        #pragma omp for reduction(+:in_circle) schedule(static)
        for (long long i = 0; i < tosses; ++i){
            lcg_state = 1664525 * lcg_state + 1013904223;
            double x = (double)(int)lcg_state * (2.0 / 4294967295.0);

            lcg_state = 1664525 * lcg_state + 1013904223;
            double y = (double)(int)lcg_state * (2.0 / 4294967295.0);

            if (x*x + y*y <= 1.0) {
                in_circle++;
            }
        }
    }
    return in_circle;
}