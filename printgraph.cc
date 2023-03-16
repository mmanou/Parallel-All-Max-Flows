/*
g++ -fopenmp -o printgraph.exe printgraph.cc
*/

#include <omp.h>
#include <climits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <string>


// #define NUM_V 4
#define NUM_V 7

int dinics_DFS(
    int curr_vertex,
    int available_flow,
    std::vector<std::vector<int>> &graph,
    std::vector<int> &levels,
    int s,
    int t
);


int main(int argc, char **argv) {

    // Read file
    std::ifstream input_file("input.csv");

    std::string line;
    std::string cell;
    std::getline(input_file, line);
    int num_v = std::stoi(line);

    std::vector<
        std::vector<int>> original_graph(num_v);

    for (int y=0; y < num_v; y++) {

        std::getline(input_file, line);
        std::stringstream line_stream(line);
        std::vector<int> curr_vector(num_v);

        for (int x=0; x < num_v; x++) {
            std::getline(line_stream, cell, ',');
            curr_vector[x] = std::stoi(cell);
        }

        original_graph[y] = curr_vector;
    }
    input_file.close();

    // Printing next graph level

    std::cout << "\nPrinting the next level of graph. num_v=" << num_v;
    for (int row=0; row < original_graph.size(); row++) {
        std::cout << "\n[";
        for (int col=0; col < original_graph[0].size(); col++) {
            std::cout << original_graph[row][col] << ", ";
        }
        std::cout << "]";
    }


    return 0;

}
