/*
g++ -fopenmp -std=c++11 -o baseline.exe baseline.cc
*/

#include <sched.h>
#include <omp.h>
#include <climits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <string>

//OMP_PLACES='{0,1}, {2,3}';

int dinics_DFS(
    int curr_vertex,
    int available_flow,
    int **graph,
    int *levels,
    int s,
    int t,
    int num_v
);

int main(int argc, char **argv) {

    double start_sequential_time = omp_get_wtime();

    // Read file
    std::ifstream input_file(argv[1]);
    std::string line;
    std::string cell;
    std::getline(input_file, line);
    const int NUM_V = std::stoi(line);

    // initialise arrays
    int **original_graph = (int **)malloc(NUM_V * sizeof(int *));
    for (int i=0; i < NUM_V; i++) {
        original_graph[i] = (int *)malloc(NUM_V * sizeof(int));
    }

    // read input file into array
    for (int y=0; y < NUM_V; y++) {
        std::getline(input_file, line);
        std::stringstream line_stream(line);
        int *curr_vector = original_graph[y];

        for (int x=0; x < NUM_V; x++) {
            std::getline(line_stream, cell, ',');
            curr_vector[x] = std::stoi(cell);
        }
    }
    input_file.close();

    const int NUM_GRAPHS = NUM_V * NUM_V - NUM_V; //number of problems to be solved in parallel

    // Create source_sink_pairs array
    int *source_sink_pairs = (int *)malloc(2*NUM_GRAPHS*sizeof(int));
    int ind=0;
    for (int i=0; i < NUM_V; i++) {
        for (int j=0; j < NUM_V; j++) {
            if (i != j) {
                source_sink_pairs[ind++] = i;
                source_sink_pairs[ind++] = j;
            }
        }
    }

    int min_max_flow = INT_MAX;

    double end_sequential_time = omp_get_wtime();

    // initialise performance data arrays
    double *perf_parallel_timings = (double *)malloc(NUM_GRAPHS*sizeof(double));
    double *perf_critical_timings = (double *)malloc(NUM_GRAPHS*sizeof(double));
    double *perf_memalloc_timings = (double *)malloc(NUM_GRAPHS*sizeof(double));
    int *perf_thread_ids = (int *)malloc(NUM_GRAPHS*sizeof(int));
    int *perf_cpu_ids = (int *)malloc(NUM_GRAPHS*sizeof(int));

    double start_parallel_time = omp_get_wtime();

    for (ind=0; ind < NUM_GRAPHS; ind++)
    {

        perf_thread_ids[ind] = omp_get_thread_num();
        perf_cpu_ids[ind] = sched_getcpu();
        perf_parallel_timings[ind] = omp_get_wtime(); //start recording time

        int s = source_sink_pairs[ind*2], t = source_sink_pairs[ind*2+1];
        std::queue<int> vertex_queue;
        int max_flow = -1;
        int curr_vertex;
        int curr_level;
        int i;

        int *levels = (int *)calloc(NUM_V, sizeof(int));

        //initialise private graph to original_graph values
        int **graph = (int **)malloc(NUM_V * sizeof(int *));
        for (int i=0; i < NUM_V; i++) {
            graph[i] = (int *)malloc(NUM_V * sizeof(int));
        }
        for (int x=0; x < NUM_V; x++) {
            for (int y=0; y < NUM_V; y++) {
                graph[x][y] = original_graph[x][y];
            }
        }

        perf_memalloc_timings[ind] = omp_get_wtime() - perf_parallel_timings[ind];

        // Remove all edges from the sink. These will sum together to become our Max Flow.
        for (i=0; i < NUM_V; i++) {
            graph[t][i] = 0;
        }

        while(max_flow == -1) {
            vertex_queue.push(s); // source vertex
            vertex_queue.push(1); // initial distance from source vertex

            // Run BFS to construct levels array, where level is min distance from s
            while (vertex_queue.size() > 0) {

                curr_vertex = vertex_queue.front();
                vertex_queue.pop();

                curr_level = vertex_queue.front();
                vertex_queue.pop();

                // Check all connected vertices for unvisited vertices
                for (i=0; i < NUM_V; i++) {

                    if (i != s
                        && i != curr_vertex
                        && graph[curr_vertex][i] > 0
                        && levels[i] == 0
                    ) {
                        // Set level
                        levels[i] = curr_level;

                        // Add non-sink vertex to queue
                        if (i != t) {
                            vertex_queue.push(i);
                            vertex_queue.push(curr_level + 1);
                        }
                    }
                }
            }

            // No route to sink. Return Max Flow.
            if (levels[t] == 0) {

                // Sum all
                max_flow = graph[t][0];
                for (i=1; i < NUM_V; i++) {
                    max_flow += graph[t][i];
                }

                perf_critical_timings[ind] = omp_get_wtime();

                // printf("Time %f Thread %d: Max Flow = %d. Waiting for critical region\n", , omp_get_thread_num(), max_flow);
                // printf("Time %f Thread %d: Enter critical region\n", omp_get_wtime(), omp_get_thread_num());
                if (max_flow < min_max_flow) {
                    min_max_flow = max_flow;
                }
                // Parallel region complete
                double curr_end_time = omp_get_wtime();
                perf_critical_timings[ind] = curr_end_time - perf_critical_timings[ind];
                perf_parallel_timings[ind] = curr_end_time - perf_parallel_timings[ind];
                //free memory
                for (i=0; i < NUM_V; i++) {
                    free(graph[i]);
                }
                free(graph);
                free(levels);

            } else {
                // Run DFS. Each step in a path moves from a lower-level to higher-level vertex
                dinics_DFS(s, INT_MAX, graph, levels, s, t, NUM_V);

                // Reset levels array
                for (i=0; i < NUM_V; i++) {
                    levels[i] = 0;
                }
            }
        }
    }

    double end_parallel_time = omp_get_wtime();

    printf("\nRESULT -- Minimum of max flows = %d --\n", min_max_flow);

    double total_sequential_time = end_sequential_time - start_sequential_time;
    double total_parallel_time = end_parallel_time - start_parallel_time;
    printf("Total Sequential time: %f\n", total_sequential_time);
    printf("Total Parallel time: %f\n", total_parallel_time);
    printf("--Parallel timings--\n");
    printf("INDEX CPU_ID THREAD_ID TOTAL_TIME MALLOC_TIME CRITICAL_TIME\n");
    // for (ind=0; ind < NUM_GRAPHS; ind++) {
    //     printf("%d %d %d %f %f %f\n",
    //         ind,
    //         perf_cpu_ids[ind],
    //         perf_thread_ids[ind],
    //         perf_parallel_timings[ind],
    //         perf_memalloc_timings[ind],
    //         perf_critical_timings[ind]
    //     );
    // }

    // free arrays
    for (int i=0; i < NUM_V; i++) {
        free(original_graph[i]);
    }
    free(original_graph);
    free(source_sink_pairs);

    free(perf_parallel_timings);
    free(perf_critical_timings);
    free(perf_memalloc_timings);
    free(perf_thread_ids);
    free(perf_cpu_ids);

    return 0;
}


int dinics_DFS(
    int curr_vertex,
    int available_flow,
    int **graph,
    int *levels,
    int s,
    int t,
    int num_v
) {
    // Base case: Sink
    if (curr_vertex == t) {
        return available_flow;
    }

    // Recursive case
    int i;
    int newly_consumed_flow;
    int total_consumed_flow = 0;
    for (i=0; i < num_v; i++) {
        // already pushed all flow
        if (available_flow == total_consumed_flow) {
            return available_flow;
        }

        // check for edge to sink or higher-level vertices
        if ((graph[curr_vertex][i] > 0) && (i == t || levels[i] > levels[curr_vertex])) {

            newly_consumed_flow = dinics_DFS(
                i,
                std::min(available_flow - total_consumed_flow, graph[curr_vertex][i]),
                graph,
                levels,
                s,
                t,
                num_v
            );

            if (newly_consumed_flow == 0) {  // No path to sink: rule out that vertex
                levels[i] = 0;
            }
            else { // Reached sink: consume flow

                graph[curr_vertex][i] -= newly_consumed_flow; // reduce forward capacity
                graph[i][curr_vertex] += newly_consumed_flow; // increase backward capacity
                total_consumed_flow += newly_consumed_flow;
            }
        }

    }
    return total_consumed_flow;
}
