#ifndef PRACTICA1_MH_PAR_ALGORITMOS_H
#define PRACTICA1_MH_PAR_ALGORITMOS_H
#include <random>
#include <cfloat>
#include "utilidades.h"
#include <chrono>

// cj es igual al cluster al que pertenece xk y la restricción es cannot-link
// o
// cj no es igual al cluster al que pertenece xk y la restrición es must-link
size_t infeasibility(size_t xi, int cj, const std::vector<int> C, const R_matrix& R) {
    size_t inf = 0;
    for (auto& v : R[xi])
        inf += (cj == C[v.first] && v.second == -1) || (cj != C[v.first] && v.second == 1);
    return inf;
}

inline size_t total_infeasibility(const std::vector<int> C, const R_matrix& R) {
    size_t total = 0;
    for (size_t xi = 0; xi < C.size(); ++xi)
        total += infeasibility(xi, C[xi], C, R);
    return total;
}

// ------------------------------GREEDY------------------------------
std::vector<int> greedy_copkm(const double_matrix& X, const R_matrix& R, std::vector<cluster>& clusters, size_t seed)
{
    std::vector<int> C(X.size(), -1); // -1 == ningún cluster, 0 == 1º cluster, 1 == 2º cluster....
    std::vector<size_t> vRSI;

    // Bajaramos los índices para recorrer X de forma aleatoria y sin repetición
    for (size_t i = 0; i < X.size(); ++i)
        vRSI.push_back(i);
    shuffle(vRSI.begin(), vRSI.end(), std::default_random_engine(seed));
    std::list<size_t> RSI(vRSI.begin(), vRSI.end());

    std::vector<int> C_ant;
    std::list<size_t> asignaciones;
    size_t inf, min, c_min;
    double d_min, d;
    for (size_t it = 0; it < 10000 && C_ant != C; ++it)
    {
        C_ant = C;

        for (auto& xi : RSI) {
            d_min = DBL_MAX;
            min = (size_t)-1;
            asignaciones.clear();

            // Calcular el incremento en infeasibility que produce la asignacion xi a cada cluster cj
            for (size_t cj = 0; cj < clusters.size(); ++cj) {
                inf = infeasibility(xi, cj, C, R);
                if (inf < min) {
                    asignaciones.clear();
                    min = inf;
                }
                if (inf <= min)
                    asignaciones.push_back(cj);
            }

            // De entre las asignaciones que produce nmenos incremento en infeasibility, seleccionar la
            // asociada al centroide más cercano a xi
            for (auto& cj : asignaciones) {
                d = clusters[cj].distancia_centroide(X[xi]);
                if (d < d_min) {
                    d_min = d;
                    c_min = cj;
                }
            }

            clusters[c_min].aniadir_punto(X[xi]);
            C[xi] = c_min;
        }

        // Actualizar el centroide promediando las instancias de su cluster asociado ci
        for (auto& ci : clusters) {
            ci.actualizar_centroide();
            ci.eliminar_puntos();
        }
    }
    
    return C;
}

void reparar_solucion(std::vector<int>& C, const R_matrix& R, size_t k)  {
    size_t x_min, inf, min;
    for (int  ci = 0; ci < k; ++ci) {
        if (count(C.begin(), C.end(), ci)==0) {
            min = (size_t)-1;
            for (size_t xi = 0; xi < R.size(); ++xi) {
                inf = infeasibility(xi, ci, C, R);
                if (inf < min && count(C.begin(), C.end(), C[xi]) > 1) {
                    min = inf;
                    x_min = xi;
                }
            }
            C[x_min] = ci;
        }
    }
}

// ------------------------------BÚSQUEDA LOCAL------------------------------
inline double fitness(const std::vector<int>& C, const double_matrix& X, const R_matrix& R, std::vector<cluster>& clusters, double lambda) {
    return desviacion_general(C, X, clusters) + lambda*total_infeasibility(C, R);
}

std::vector<int> busqueda_trayectorias_simples(const double_matrix& X, const R_matrix& R, std::vector<cluster>& clusters, double lambda)
{
    std::vector<int> C, S;
    size_t k = clusters.size(), n = X.size();

    // Generación de la solución inicial
    while (empty_clusters(C, k) > 0) {
        C.clear();
        for (size_t i = 0; i < n; ++i)
            C.push_back(rand()%k);
    }

    double f_actual = fitness(C, X, R, clusters, lambda), f_vecino;
    bool hay_mejora = true;
    for (size_t it = 0; it < 10000 && hay_mejora; ++it) {
        hay_mejora = false;

        // Exploración aleatoria del entorno y selección del primer mejor vecino
        for (size_t i = rand()%n, ci = 0; ci < n && !hay_mejora; i = (i+1)%n, ++ci) {
            for (size_t l = rand()%k, cl = 0; cl < k && !hay_mejora; l = (l+1)%k, ++cl) {
                S = C;
                S[i] = l;
                f_vecino = fitness(S, X, R, clusters, lambda);
                if (f_vecino <= f_actual && empty_clusters(S, k) == 0) {
                    f_actual = f_vecino;
                    hay_mejora = true;
                    C = S;
                }
            }
        }
    }
    std::cout << "Agregado: " << f_actual << std::endl;
    return C;
}

#endif //PRACTICA1_MH_PAR_ALGORITMOS_H