#include <iostream>
#include "algoritmos.h"
#include <chrono>
#define ZOO_10 486
#define ZOO_20 912
#define GLASS_10 2167
#define GLASS_20 4096
#define BUPA_10 5632
#define BUPA_20 10824

// zoo seeds: 1618392344, 1618393097, 1618393257, 1618393407, 1618393577
// glass seeds: 1618395394, 1618397214, 1618397437, 1618397651, 1618398759
// bupa seeds: 1618398086, 1618398519, 1618398924, 1618399070, 1618399285
int main()
{
    // Inicializar semillas y número de clusters
    size_t seed = time(NULL), k = 7;
    srand(seed);

    // Leer datos
    R_matrix R = matriz_a_lista(archivo_a_matriz<int>(std::ifstream("data/zoo_set_const_20.const")));
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

    std::cout << "-Búsqueda por trayectorias simples-" << std::endl;
    // Ejecución de búsqueda local
    double lambda = distancia_maxima(X)*10/ZOO_20;
    start_time = std::chrono::system_clock::now();
    C = busqueda_trayectorias_simples(X, R, clusters, lambda);
    end_time = std::chrono::system_clock::now();

    std::cout << "Infactibilidad: " << total_infeasibility(C, R) << std::endl;
    std::cout << "Desviación general: " << desviacion_general(C, X, clusters) << std::endl;
    std::cout << "Clusters vacíos: " << empty_clusters(C, k) << std::endl;
    std::cout << "Tiempo transcurrido: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms " << std::endl;
    std::cout << "Semilla: " << seed << std::endl << std::endl;

    return 0;
}
