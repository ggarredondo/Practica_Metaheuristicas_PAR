#include <iostream>
#include "algoritmos.h"

//      seeds
// zoo: 2, 7
int main() {
    size_t seed = time(NULL), k = 16;
    srand(seed);
    R_matrix R = matriz_a_lista(archivo_a_matriz<int>(std::ifstream("data/zoo_set_const_10.const")));
    double_matrix X = archivo_a_matriz<double>(std::ifstream("data/zoo_set.dat"));
    std::vector<cluster> clusters;
    for (size_t i = 0; i < k; ++i)
        clusters.push_back(cluster(seed, X[0].size()));
    std::vector<int> C = greedy_copkm(X, R, clusters, seed);
    reparar_solucion(C, R, k);

    std::cout << "Infactibilidad: " << total_infeasibility(C, R) << std::endl;
    std::cout << "Clusters vacíos: " << empty_clusters(C, k) << std::endl;
    std::cout << seed << std::endl;

    return 0;
}
