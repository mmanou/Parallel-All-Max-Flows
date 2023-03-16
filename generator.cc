/*
g++ -o generator.exe generator.cc
generator.exe ./_filename_.csv
*/

#include <cstdlib>
#include <ctime>
#include <string>
#include <fstream>
#include <iostream>

void generate_graph_file(
        int n,
        int edge_probability,
        int max_edge_weight,
        std::string filename);


int main(int argc, char **argv) {
    srand(time(0));

    std::string filename = argv[1];

    std::cout << "Input n, then p (where prob. edge = 1/p), then max_edge_weight\n";

    int n, edge_probability, max_edge_weight;
    std::cin >> n >> edge_probability >> max_edge_weight;

    generate_graph_file(n, edge_probability, max_edge_weight - 1, filename);

    return 0;
}


void generate_graph_file(int n, int edge_probability, int max_edge_weight, std::string filename) {

    int num_edges = 0;

    std::ofstream outfile(filename);

    outfile << n << "\n";

    for (int y=0; y < n; y++) {
        for (int x=0; x < n; x++) {
            if (x == y) {
                // No arrow to self
                outfile << 0;
            } else if (
                    x == y+1 || // guaranteed arrow
                    (y == n-1 && x == 0) || // guaranteed arrow on last line
                    (rand() % edge_probability == 0)) { // random arrow
                num_edges++;
                outfile << ((rand() % max_edge_weight) + 1);
            } else  {
                // random non-arrow
                outfile << 0;
            }

            if (x < n-1) {
                outfile << ",";
            }
        }
        outfile << "\n";
    }

    std::cout << num_edges << "\n";
}
