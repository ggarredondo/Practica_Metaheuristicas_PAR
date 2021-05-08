#include <iostream>
#include "algoritmos.h"
#include <chrono>
#include <string>
#define ZOO_10 486
#define ZOO_20 912
#define GLASS_10 2167
#define GLASS_20 4096
#define BUPA_10 5632
#define BUPA_20 10824

bool preparar_datos(const std::string& set, size_t res, std::string& X_file, std::string& R_file, size_t& k) {
    bool exito = false;
    if (set.compare("bupa") == 0 && (res == 10 || res == 20)) {
        k = 16;
        exito = true;
    }
    else if ((set.compare("zoo") == 0 || set.compare("glass") == 0) && (res == 10 || res == 20)) {
        k = 7;
        exito = true;
    }
    std::string base = std::string("data/").append(set);
    X_file = base;
    X_file.append("_set.dat");
    R_file = base;
    R_file.append("_set_const_");
    R_file.append(std::to_string(res));
    R_file.append(".const");
    return exito;
}

void error() {
    std::cout << "Error en el formato." << std::endl;
    std::cout << "formato: ./p1_par <nombre_set> <porcentaje_rest> <semilla>" << std::endl;
    std::cout << "\n-<nombre_set> debe estar en minúsculas." << std::endl;
    std::cout << "-<porcentaje_rest> debe ser 10 o 20." << std::endl;
}

// formato: ./p1_par <nombre_set> <porcentaje_rest> <semilla>
// zoo seeds: 1618392344, 1618393097, 1618393257, 1618393407, 1618393577
// glass seeds: 1618395394, 1618397214, 1618397437, 1618397651, 1618398759
// bupa seeds: 1618398086, 1618398519, 1618398924, 1618399070, 1618399285
int main(int argc, char *argv[])
{
    if (argc != 4) {
        error();
        return -1;
    }
    // inicializar parámetros y leer argumentos
    std::string set = argv[1], X_file, R_file;
    size_t k, res = std::stoi(argv[2]), seed = std::stoi(argv[3]);
    srand(seed);
    if (preparar_datos(set, res, X_file, R_file, k))
        std::cout << "Formato correcto.\n" << std::endl;
    else {
        error();
        return -1;
    }

    // Leer datos
    int_matrix m = archivo_a_matriz<int>(std::ifstream(R_file));
    R_matrix R = matriz_a_lista(m);
    R_list Rlista = matriz_a_Rlista(m);
    double_matrix X = archivo_a_matriz<double>(std::ifstream(X_file));
    std::vector<cluster> clusters;
    double lambda = distancia_maxima(X)/Rlista.size();

    // Inicializar clusters
    for (size_t i = 0; i < k; ++i)
        clusters.push_back(cluster(seed, X[0].size()));

    // Ejecución de greedy
    std::chrono::time_point<std::chrono::system_clock> start_time = std::chrono::system_clock::now();
    std::vector<int> C = greedy_copkm(X, R, clusters, seed);
    reparar_solucion(C, R, k);
    std::chrono::time_point<std::chrono::system_clock> end_time = std::chrono::system_clock::now();

    std::cout << "-Greedy-\nAgregado: " << fitness(C, X, Rlista, clusters, lambda) << std::endl;
    std::cout << "Infactibilidad: " << total_infeasibility(C, Rlista) << std::endl;
    std::cout << "Desviación general: " << desviacion_general(C, X, clusters) << std::endl;
    std::cout << "Tiempo transcurrido: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms " << std::endl;
    std::cout << "Semilla: " << seed << std::endl << std::endl;

    // Ejecución de búsqueda local
    start_time = std::chrono::system_clock::now();
    C = busqueda_local(X, Rlista, clusters, lambda, seed);
    end_time = std::chrono::system_clock::now();

    std::cout << "-Búsqueda local-\nAgregado: " << fitness(C, X, Rlista, clusters, lambda) << std::endl;
    std::cout << "Infactibilidad: " << total_infeasibility(C, Rlista) << std::endl;
    std::cout << "Desviación general: " << desviacion_general(C, X, clusters) << std::endl;
    std::cout << "Tiempo transcurrido: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms " << std::endl;
    std::cout << "Semilla: " << seed << std::endl << std::endl;

    // Ejecución de algoritmo genético generacional con cruce informe
    start_time = std::chrono::system_clock::now();
    C = AGG_UN(X, Rlista, clusters, lambda, seed);
    reparar_solucion(C, R, k);
    end_time = std::chrono::system_clock::now();

    std::cout << "-AGG UN-\nAgregado: " << fitness(C, X, Rlista, clusters, lambda) << std::endl;
    std::cout << "Infactibilidad: " << total_infeasibility(C, Rlista) << std::endl;
    std::cout << "Desviación general: " << desviacion_general(C, X, clusters) << std::endl;
    std::cout << "Tiempo transcurrido: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms " << std::endl;
    std::cout << "Semilla: " << seed << std::endl << std::endl;

    return 0;
}
