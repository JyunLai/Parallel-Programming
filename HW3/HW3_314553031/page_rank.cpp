#include "page_rank.h"

#include <cmath>
#include <cstdlib>
#include <omp.h>

#include "../common/graph.h"

// page_rank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is num_nodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void page_rank(Graph g, double *solution, double damping, double convergence)
{

    // initialize vertex weights to uniform probability. Double
    // precision scores are used to avoid underflow for large graphs

    int nnodes = num_nodes(g);
    double equal_prob = 1.0 / nnodes;
    double *new_score = new double[nnodes];
    #pragma omp parallel for
    for (int i = 0; i < nnodes; ++i)
    {
        solution[i] = equal_prob;
    }

    bool cover = false;
    const double none_prob = (1.0 - damping) / nnodes;
    while (!cover) {
      double diff = 0.0;

      double contribute = 0.0;
      #pragma omp parallel for reduction(+:contribute)
      for (int i = 0; i < nnodes; i++) {
        if (outgoing_size(g, i) == 0) {
          contribute += damping * solution[i] / nnodes;
        }
      }

      #pragma omp parallel for
      for (int i = 0; i < nnodes; i++) {
        double sum = 0.0;
        const Vertex *start = incoming_begin(g, i);
        const Vertex *end = incoming_end(g, i);
        for (const Vertex *ptr = start; ptr != end; ptr++) {
          int tmp = *ptr;
          int degree = outgoing_size(g, tmp);
          if (degree > 0) {
            sum += solution[tmp] / (double)degree;
          }
        }
        new_score[i] = (damping * sum) + none_prob;
        new_score[i] += contribute;
      }

      #pragma omp parallel for reduction(+:diff)
      for (int i = 0; i < nnodes; i++) {
        diff += std::abs(new_score[i] - solution[i]);
        solution[i] = new_score[i];
      }

      cover = (diff < convergence);
    }
    delete[] new_score;
}
