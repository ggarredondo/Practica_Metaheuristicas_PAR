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
    mostrar_resultados("Greedy", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de búsqueda local
    start_time = std::chrono::system_clock::now();
    C = busqueda_local(X, Rlista, clusters, lambda, seed);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("Búsqueda local", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de algoritmo genético generacional con cruce uniforme
    start_time = std::chrono::system_clock::now();
    C = AGG_UN(X, Rlista, clusters, lambda, seed);
    reparar_solucion(C, R, k);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("AGG_UN", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de algoritmo genético generacional con cruce por segmento fijo
    start_time = std::chrono::system_clock::now();
    C = AGG_SF(X, Rlista, clusters, lambda, seed);
    reparar_solucion(C, R, k);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("AGG_SF", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de algoritmo genético estacionario con cruce uniforme
    start_time = std::chrono::system_clock::now();
    C = AGE_UN(X, Rlista, clusters, lambda, seed);
    reparar_solucion(C, R, k);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("AGE_UN", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de algoritmo genético estacionario con cruce por segmento fijo
    start_time = std::chrono::system_clock::now();
    C = AGE_SF(X, Rlista, clusters, lambda, seed);
    reparar_solucion(C, R, k);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("AGE_SF", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de algoritmo memético con probabilidad de explotación 1
    start_time = std::chrono::system_clock::now();
    C = AM(X, Rlista, clusters, lambda, seed, 1.0f, false);
    reparar_solucion(C, R, k);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("AM_(10,1.0)", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de algoritmo memético con probabilidad de explotación 0.1
    start_time = std::chrono::system_clock::now();
    C = AM(X, Rlista, clusters, lambda, seed, 0.1f, false);
    reparar_solucion(C, R, k);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("AM_(10,0.1)", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    // Ejecución de algoritmo memético con explotación a los 0.1*N mejores
    start_time = std::chrono::system_clock::now();
    C = AM(X, Rlista, clusters, lambda, seed, 0.1f, true);
    reparar_solucion(C, R, k);
    end_time = std::chrono::system_clock::now();
    mostrar_resultados("AM_(10,0.1mej)", fitness(C, X, Rlista, clusters, lambda), total_infeasibility(C, Rlista), desviacion_general(C, X, clusters), std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count());

    std::cout << "Semilla: " << seed << std::endl;

    return 0;
}
