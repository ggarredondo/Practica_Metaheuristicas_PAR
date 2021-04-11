#include <iostream>
#include "algoritmos.h"
#include <chrono>
#define ZOO_10 486
#define ZOO_20 912
#define GLASS_10 2167
#define GLASS_20 4096
#define BUPA_10 5632
#define BUPA_20 10824

// zoo seeds to test: 1617873919, 1617873954
// glass seeds to test: 1617874026, 1617874083
// bupa seeds to test: 1617872058, 1617873434
int main()
{
    // Inicializar semillas y número de clusters
    size_t seed = time(NULL), k = 7;
    srand(seed);

    // Leer datos
    R_matrix R = matriz_a_lista(archivo_a_matriz<int>(std::ifstream("data/zoo_set_const_10.const")));
    double_matrix X = archivo_a_matriz<double>(std::ifstream("data/zoo_set.dat"));
    std::vector<cluster> clusters;

    // Inicializar clusters
    for (size_t i = 0; i < k; ++i)
        clusters.push_back(cluster(seed, X[0].size()));

    // Ejecución de greedy
    std::chrono::time_point<std::chrono::system_clock> start_time = std::chrono::system_clock::now();
    std::vector<int> C = greedy_copkm(X, R, clusters, seed);
    reparar_solucion(C, R, k);
    std::chrono::time_point<std::chrono::system_clock> end_time = std::chrono::system_clock::now();

    std::cout << "-Greedy-\nInfactibilidad: " << total_infeasibility(C, R) << std::endl;
    std::cout << "Desviación general: " << desviacion_general(C, X, clusters) << std::endl;
    std::cout << "Clusters vacíos: " << empty_clusters(C, k) << std::endl;
    std::cout << "Tiempo transcurrido: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms " << std::endl;
    std::cout << "Semilla: " << seed << std::endl << std::endl;

    // Ejecución de búsqueda local
    double lambda = distancia_maxima(X)/ZOO_10;
    start_time = std::chrono::system_clock::now();
    C = busqueda_trayectorias_simples(X, R, clusters, lambda);
    end_time = std::chrono::system_clock::now();

    std::cout << "-Búsqueda por trayectorias simples-\nInfactibilidad: " << total_infeasibility(C, R) << std::endl;
    std::cout << "Desviación general: " << desviacion_general(C, X, clusters) << std::endl;
    std::cout << "Clusters vacíos: " << empty_clusters(C, k) << std::endl;
    std::cout << "Tiempo transcurrido: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms " << std::endl;
    std::cout << "Semilla: " << seed << std::endl << std::endl;

    return 0;
}
