#include "bfs.h"
#include <cstring>
#include <cstdlib>
#include <omp.h>

#include "../common/graph.h"

#ifdef VERBOSE
#include "../common/CycleTimer.h"
#include <stdio.h>
#endif // VERBOSE

constexpr int ROOT_NODE_ID = 0;
constexpr int NOT_VISITED_MARKER = -1;

void vertex_set_clear(VertexSet *list)
{
    list->count = 0;
}

void vertex_set_init(VertexSet *list, int count)
{
    list->max_vertices = count;
    list->vertices = new int[list->max_vertices];
    vertex_set_clear(list);
}

void vertex_set_destroy(VertexSet *list)
{
    delete[] list->vertices;
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(Graph g, VertexSet *frontier, VertexSet *new_frontier, int *distances, VertexSet* local_frontiers)
{
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < frontier->count; i++)
    {
        int thread_id = omp_get_thread_num();
        VertexSet* my_frontier = &local_frontiers[thread_id];

        int node = frontier->vertices[i];
        int new_distance = distances[node] + 1;
        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1) ? g->num_edges : g->outgoing_starts[node + 1];

        // attempt to add all neighbors to the new frontier
        for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
        {
            int outgoing = g->outgoing_edges[neighbor];

            if (__sync_bool_compare_and_swap(&distances[outgoing], NOT_VISITED_MARKER, new_distance)) {
                int index = my_frontier->count++;
                my_frontier->vertices[index] = outgoing;
            }
        }
    }
}

void bottom_up_step(Graph g, VertexSet *new_frontier, int *distances, int current_distance, VertexSet* local_frontiers)
{
    int total_nodes = num_nodes(g);
    int new_distance = current_distance + 1;

    #pragma omp parallel for schedule(dynamic, 512)
    for (int i = 0; i < total_nodes; i++) {
        int thread_id = omp_get_thread_num();
        VertexSet* my_frontier = &local_frontiers[thread_id];

        if (distances[i] == NOT_VISITED_MARKER) {
            const Vertex* start = incoming_begin(g, i);
            const Vertex* end = incoming_end(g, i);

            for (const Vertex *ptr = start; ptr != end; ptr++) {
                int j = *ptr;
                if (distances[j] == current_distance) {
                    if (__sync_bool_compare_and_swap(&distances[i], NOT_VISITED_MARKER, new_distance)) {
                        int index = my_frontier->count++;
                        my_frontier->vertices[index] = i;
                    }
                    break;
                }
            }
        }
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution *sol)
{

    VertexSet list1;
    VertexSet list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    VertexSet *frontier = &list1;
    VertexSet *new_frontier = &list2;

    int max_threads = omp_get_max_threads();
    VertexSet* local_frontiers = new VertexSet[max_threads];
    for (int i = 0; i < max_threads; i++) {
        vertex_set_init(&local_frontiers[i], graph->num_nodes);
    }

    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0)
    {

#ifdef VERBOSE
        double start_time = CycleTimer::current_seconds();
#endif

        vertex_set_clear(new_frontier);
        for (int i = 0; i < max_threads; i++) {
            vertex_set_clear(&local_frontiers[i]);
        }

        top_down_step(graph, frontier, new_frontier, sol->distances, local_frontiers);

        int* offset = new int[max_threads + 1];
        offset[0] = 0;
        for (int i = 0; i < max_threads; i++) {
            offset[i+1] = offset[i] + local_frontiers[i].count;
        }

        int global_count = offset[max_threads];
        new_frontier->count = global_count;
        #pragma omp parallel for
        for (int i = 0; i < max_threads; i++) {
            int my_offset = offset[i];
            VertexSet* my_frontier = &local_frontiers[i];
            std::memcpy(new_frontier->vertices + my_offset, my_frontier->vertices, my_frontier->count * sizeof(int));
        }
        delete[] offset;

#ifdef VERBOSE
        double end_time = CycleTimer::current_seconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        VertexSet *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }

    for (int i = 0; i < max_threads; i++) {
        vertex_set_destroy(&local_frontiers[i]);
    }
    delete[] local_frontiers;

    // free memory
    vertex_set_destroy(&list1);
    vertex_set_destroy(&list2);
}

void bfs_bottom_up(Graph graph, solution *sol)
{
    VertexSet list1;
    VertexSet list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    VertexSet *frontier = &list1;
    VertexSet *new_frontier = &list2;

    int max_threads = omp_get_max_threads();
    VertexSet* local_frontiers = new VertexSet[max_threads];
    for (int i = 0; i< max_threads; i++) {
        vertex_set_init(&local_frontiers[i], graph->num_nodes);
    }

    int distance = 0;
    #pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::current_seconds();
#endif

        vertex_set_clear(new_frontier);
        for (int i = 0; i < max_threads; i++) {
            vertex_set_clear(&local_frontiers[i]);
        }

        bottom_up_step(graph, new_frontier, sol->distances, distance, local_frontiers);

        int* offset = new int[max_threads + 1];
        offset[0] = 0;
        for (int i = 0; i < max_threads; i++) {
            offset[i+1] = offset[i] + local_frontiers[i].count;
        }

        int global_count = offset[max_threads];
        new_frontier->count = global_count;
        #pragma omp parallel for
        for (int i = 0; i < max_threads; i++) {
            int my_offset = offset[i];
            VertexSet* my_frontier = &local_frontiers[i];
            std::memcpy(new_frontier->vertices + my_offset, my_frontier->vertices, my_frontier->count * sizeof(int));
        }
        delete[] offset;

#ifdef VERBOSE
        double end_time = CycleTimer::current_seconds();
        printf("frontier=%-10d %.4f sec\n", new_frontier->count, end_time - start_time);
#endif

        VertexSet *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
        distance++;
    }
    for (int i = 0; i < max_threads; i++) {
        vertex_set_destroy(&local_frontiers[i]);
    }
    delete[] local_frontiers;

    vertex_set_destroy(&list1);
    vertex_set_destroy(&list2);
}

void bfs_hybrid(Graph graph, solution *sol)
{
    const int top_down_threshold = graph->num_nodes / 20;

    VertexSet list1;
    VertexSet list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    VertexSet *frontier = &list1;
    VertexSet *new_frontier = &list2;

    int max_threads = omp_get_max_threads();
    VertexSet* local_frontiers = new VertexSet[max_threads];
    for (int i = 0; i < max_threads; i++) {
        vertex_set_init(&local_frontiers[i], graph->num_nodes);
    }

    int distance = 0;
    #pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::current_seconds();
#endif

        vertex_set_clear(new_frontier);
        for (int i = 0; i < max_threads; i++) {
            vertex_set_clear(&local_frontiers[i]);
        }

        if (frontier->count > top_down_threshold) {

#ifdef VERBOSE
            printf("(Botton-Up)");
#endif

            bottom_up_step(graph, new_frontier, sol->distances, distance, local_frontiers);
        }
        else {

#ifdef VERBOSE
            printf("(Top-Down)");
#endif

            top_down_step(graph, frontier, new_frontier, sol->distances, local_frontiers);
        }

        int* offset = new int[max_threads + 1];
        offset[0] = 0;
        for (int i = 0; i < max_threads; i++) {
            offset[i+1] = offset[i] + local_frontiers[i].count;
        }

        int global_count = offset[max_threads];
        new_frontier->count = global_count;
        #pragma omp parallel for
        for (int i = 0; i < max_threads; i++) {
            int my_offset = offset[i];
            VertexSet* my_frontier = &local_frontiers[i];
            std::memcpy(new_frontier->vertices + my_offset, my_frontier->vertices, my_frontier->count * sizeof(int));
        }
        delete[] offset;

#ifdef VERBOSE
        double end_time = CycleTimer::current_seconds();
        printf("frontier=%-10d %.4f sec\n", new_frontier->count, end_time - start_time);
#endif

        VertexSet *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
        distance++;
    }
    for (int i = 0; i < max_threads; i++) {
        vertex_set_destroy(&local_frontiers[i]);
    }
    delete[] local_frontiers;

    vertex_set_destroy(&list1);
    vertex_set_destroy(&list2);
}
