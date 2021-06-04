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

void results_to_file(const std::string& set, size_t res) {
    std::string output = "salida/";
    output.append(set);
    output.append(std::to_string(res));
    output.append(".csv");
    std::ofstream file(output, std::ios::out | std::ios::trunc);
    std::vector<size_t> zoo_seeds = {1618392344, 1618393097, 1618393257, 1618393407, 1618393577};
    std::vector<size_t> glass_seeds = {1618395394, 1618397214, 1618397437, 1618397651, 1618398759};
    std::vector<size_t> bupa_seeds = {1618398086, 1618398519, 1618398924, 1618399070, 1618399285};
    std::vector<size_t> chosen;
    if (set.compare("zoo") == 0)
        chosen = zoo_seeds;
    else if (set.compare("glass") == 0)
        chosen = glass_seeds;
    else if (set.compare("bupa") == 0)
        chosen = bupa_seeds;
    std::string X_file, R_file;
    size_t k;
    // Leer datos
    preparar_datos(set, res, X_file, R_file, k);
    int_matrix m = archivo_a_matriz<int>(std::ifstream(R_file));
    R_matrix R = matriz_a_lista(m);
    R_list Rlista = matriz_a_Rlista(m);
    double_matrix X = archivo_a_matriz<double>(std::ifstream(X_file));
    size_t n = X.size(), n_seeds = chosen.size();
    std::vector<cluster> clusters;
    double lambda = distancia_maxima(X)/Rlista.size();

    std::chrono::time_point<std::chrono::system_clock> start_time;
    std::chrono::time_point<std::chrono::system_clock> end_time;
    std::vector<int> C;

    file << "Greedy" << std::endl;
    for (size_t s = 0; s < n_seeds; ++s) {
        // Ejecución de Greedy
        srand(chosen[s]);
        // Inicializar clusters
        clusters.clear();
        for (size_t i = 0; i < k; ++i)
            clusters.push_back(cluster(X[0].size()));
        start_time = std::chrono::system_clock::now();
        C = greedy_copkm(X, R, clusters, chosen[s]);
        reparar_solucion(C, R, k);
        end_time = std::chrono::system_clock::now();
        file << total_infeasibility(C, Rlista) << "," << desviacion_general(C, X, clusters) << "," << fitness(C, X, Rlista, clusters, lambda) << "," << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()*0.001 << std::endl;
    }

    file << "BL" << std::endl;
    for (size_t s = 0; s < n_seeds; ++s) {
        // Ejecución de búsqueda local
        srand(chosen[s]);
        start_time = std::chrono::system_clock::now();
        C = busqueda_local(generar_solucion_aleatoria(n, k), X, Rlista, clusters, lambda, max_evaluaciones, chosen[s]);
        end_time = std::chrono::system_clock::now();
        file << total_infeasibility(C, Rlista) << "," << desviacion_general(C, X, clusters) << "," << fitness(C, X, Rlista, clusters, lambda) << "," << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()*0.001 << std::endl;
    }

    file << "ES" << std::endl;
    for (size_t s = 0; s < n_seeds; ++s) {
        // Ejecución de enfriamiento simulado
        srand(chosen[s]);
        start_time = std::chrono::system_clock::now();
        C = enfriamiento_simulado(generar_solucion_aleatoria(n, k), X, Rlista, clusters, lambda, max_evaluaciones, chosen[s]);
        end_time = std::chrono::system_clock::now();
        file << total_infeasibility(C, Rlista) << "," << desviacion_general(C, X, clusters) << "," << fitness(C, X, Rlista, clusters, lambda) << "," << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()*0.001 << std::endl;
    }

    file << "BMB" << std::endl;
    for (size_t s = 0; s < n_seeds; ++s) {
        // Ejecución de búsqueda multiarranque básica
        srand(chosen[s]);
        start_time = std::chrono::system_clock::now();
        C = BMB(X, Rlista, clusters, lambda, chosen[s]);
        end_time = std::chrono::system_clock::now();
        file << total_infeasibility(C, Rlista) << "," << desviacion_general(C, X, clusters) << "," << fitness(C, X, Rlista, clusters, lambda) << "," << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()*0.001 << std::endl;
    }

    file << "ILS" << std::endl;
    for (size_t s = 0; s < n_seeds; ++s) {
        // Ejecución de búsqueda local iterativa
        srand(chosen[s]);
        start_time = std::chrono::system_clock::now();
        C = ILS(X, Rlista, clusters, lambda, chosen[s]);
        end_time = std::chrono::system_clock::now();
        file << total_infeasibility(C, Rlista) << "," << desviacion_general(C, X, clusters) << "," << fitness(C, X, Rlista, clusters, lambda) << "," << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()*0.001 << std::endl;
    }

    file << "ILS_ES" << std::endl;
    for (size_t s = 0; s < n_seeds; ++s) {
        // Ejecución de ILS-ES
        srand(chosen[s]);
        start_time = std::chrono::system_clock::now();
        C = ILS_ES(X, Rlista, clusters, lambda, chosen[s]);
        end_time = std::chrono::system_clock::now();
        file << total_infeasibility(C, Rlista) << "," << desviacion_general(C, X, clusters) << "," << fitness(C, X, Rlista, clusters, lambda) << "," << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()*0.001 << std::endl;
    }
}

void error() {
    std::cout << "Error en el formato." << std::endl;
    std::cout << "formato: ./p1_par <nombre_set> <porcentaje_rest> <semilla>" << std::endl;
    std::cout << "\n-<nombre_set> debe estar en minúsculas." << std::endl;
    std::cout << "-<porcentaje_rest> debe ser 10 o 20." << std::endl;
}

void mostrar_resultados(const std::string& algoritmo, double fitness, size_t infactibilidad, double desviacion, double tiempo) {
    std::cout << "-" << algoritmo << "-" << std::endl;
    std::cout << "Infactibilidad: " << infactibilidad << std::endl;
    std::cout << "Desviación general: " << desviacion << std::endl;
    std::cout << "Agregado: " << fitness << std::endl;
    std::cout << "Tiempo transcurrido: " << tiempo << " ms " << std::endl << std::endl;
}

// formato: ./p1_par <nombre_set> <porcentaje_rest> <semilla>
// zoo seeds: 1618392344, 1618393097, 1618393257, 1618393407, 1618393577
// glass seeds: 1618395394, 1618397214, 1618397437, 1618397651, 1618398759
// bupa seeds: 1618398086, 1618398519, 1618398924, 1618399070, 1618399285
int main(int argc, char *argv[])
{
    //results_to_file("bupa", 20);
    //std::cout << "finished" << std::endl;
    //return 0;
    if (argc != 4) {
        error();
        return -1;
    }
    // inicializar parámetros y leer argumentos
    std::string set = argv[1], X_file, R_file;
    size_t k, res = std::stoi(argv[2]), seed = std::stoi(argv[3]);
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
    size_t n = X.size();
    std::vector<cluster> clusters;
    double lambda = distancia_maxima(X)/Rlista.size();

    srand(seed);
    // Inicializar clusters
    for (size_t i = 0; i < k; ++i)
        clusters.push_back(cluster(X[0].size()));

    // Ejecución de greedy
    std::chrono::time_point<std::chrono::system_clock> start_time = std::chrono::system_clock::now();
    std::vector<int> C = greedy_copkm(X, R, clusters, seed);
    reparar_solucion(C, R, k);
    std::chrono::time_point<std::chrono::system_clock> end_time = std::chrono::system_clock::now();
    mostrar_resultados("Greedy", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de búsqueda local
    srand(seed);
    start_time = std::chrono::system_clock::now();
    C = busqueda_local(generar_solucion_aleatoria(n, k), X, Rlista, clusters, lambda, max_evaluaciones, seed);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("Búsqueda local", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de enfriamiento simulado
    srand(seed);
    start_time = std::chrono::system_clock::now();
    C = enfriamiento_simulado(generar_solucion_aleatoria(n, k), X, Rlista, clusters, lambda, max_evaluaciones, seed);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("Enfriamiento simulado", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de búsqueda multiarranque básica
    srand(seed);
    start_time = std::chrono::system_clock::now();
    C = BMB(X, Rlista, clusters, lambda, seed);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("BMB", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de búsqueda local iterativa
    srand(seed);
    start_time = std::chrono::system_clock::now();
    C = ILS(X, Rlista, clusters, lambda, seed);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("ILS", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de ILS-ES
    srand(seed);
    start_time = std::chrono::system_clock::now();
    C = ILS_ES(X, Rlista, clusters, lambda, seed);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("ILS_ES", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());


    std::cout << "Semilla: " << seed << std::endl;

    return 0;
}
