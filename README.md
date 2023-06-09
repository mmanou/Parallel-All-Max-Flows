# Parallel All Max Flows

This program performs high-performance computation of the following problem:

* A ___flow network___ is a weighted directed graph with one source vertex $s$ and one sink vertex $t$.
* The ___maximum flow___ of a flow network is the ___minimum total weight of edges___ which needs to be removed in order for there to be ___no route___ from vertex $s$ to vertex $t$.
* Find the ___minimum___ of the maximum flows for the set of all flow networks that can be defined given a ___weighted directed graph with no self-loops___. That is, for all possible pairs $(s, t)$ such that $s \neq t$. All flow networks in this problem are expected to have at least one valid route from $s$ to $t$. Therefore, each flow network has a maximum flow $> 0$.

## Algorithm

The algorithm is optimised for "dense graph" inputs, where $E > V^2$. With this in mind, the following decisions were made regarding algorithm design:
- _Dinitz's algorithm_ (also known as _Dinic's algorithm_) was chosen for the sequential part of the algorithm, as it has a time complexity of $O(V^2 E)$. The other candidate for a maximum flow sequential algorithm was the _Edmonds-Karp algorithm_, with a time complexity of $O(V E^2)$.
- Graphs were represented using adjacency matrices rather than adjacency lists, as reading, creating, and deleting an edge can be performed in constant time.
- A Depth-First Search algorithm was selected for the blocking flow search part of Dinitz's Algorithm, as this has a time complexity of $O(V)$. The other candidate was an Advance and Retreat method, with a time complexity of $O(E)$.

## Testing

Tests were run on a single node on the ‘physical’ partition of the Spartan High Performance Computing system, operated by Research Computing Services at The University of Melbourne (https://dashboard.hpc.unimelb.edu.au/).
The tested node is comprised of 72 cores split evenly among 4 NUMA nodes, exclusively comprised of Intel® Xeon® Gold 6254 Processors (24.75M Cache, 3.10 GHz).

## Files

|  Filename            |  Description                                   |
|---------------------:|:-----------------------------------------------|
| solution.cc          | Optimised solution.                            |
| solution_critical.cc | Previous version. Updates best result in critical section, rather than array in parallel. |
| baseline.cc          | Runs algorithm on single core, no parallelism. |
| generator.cc         | Generates an arbitrary graph for use as input. |
| printgraph.cc        | Displays an input graph as text.               |

## Running Instructions
* Generate an input graph: 
```
g++ -o generator.exe generator.cc
./generator.exe ./<filename>.csv
```

* (Optional) Print the graph matrix to terminal:
```
g++ -o printgraph.exe printgraph.cc
./printgraph.exe ./<filename>.csv
```

* Run the solver:
```
g++ -fopenmp -std=c++11 -o solution.exe solution.cc
./solution.exe ./<filename>.csv
```
