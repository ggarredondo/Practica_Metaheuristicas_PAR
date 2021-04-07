#include <iostream>
#include "algoritmos.h"

int main() {
    size_t seed = 1, k = 16;
    srand(seed);
    R_matrix R = matriz_a_lista(archivo_a_matriz<int>(std::ifstream("data/bupa_set_const_10.const")));
    double_matrix X = archivo_a_matriz<double>(std::ifstream("data/bupa_set.dat"));
    std::vector<cluster> clusters;
    for (size_t i = 0; i < k; ++i)
        clusters.push_back(cluster(seed, X[0].size()));
    std::vector<int> C = greedy_copkm(X, R, clusters, seed);

    for (auto &v : C)
        std::cout << v << std::endl;

    return 0;
}
