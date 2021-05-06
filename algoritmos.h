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

inline size_t total_infeasibility(const std::vector<int>& C, const R_list& R) {
    size_t total = 0;
    for (auto& r : R)
        total += (C[r.i] == C[r.j] && r.r == -1) ||  (C[r.i] != C[r.j] && r.r == 1);
    return total;
}

inline double fitness(const std::vector<int>& C, const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda) {
    return desviacion_general(C, X, clusters) + lambda*total_infeasibility(C, R);
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
std::vector<int> busqueda_local(const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t seed)
{
    std::vector<int> C, S;
    std::vector<std::pair<size_t, size_t> > vecindario;
    size_t k = clusters.size(), n = X.size(), aux;

    // Generación de la solución inicial
    while (empty_clusters(C, k) > 0) {
        C.clear();
        for (size_t i = 0; i < n; ++i)
            C.push_back(rand()%k);
    }

    double f_actual = fitness(C, X, R, clusters, lambda), f_vecino;
    bool hay_mejora = true;
    for (size_t ev = 0; ev < 100000 && hay_mejora; ++ev) {
        hay_mejora = false;

        // Generación del entorno
        vecindario.clear();
        for (size_t i = 0; i < n; ++i) {
            for (size_t l = 0; l < k; ++l) {
                if (C[i] != l) {
                    aux = C[i];
                    C[i] = l;
                    if (empty_clusters(C, k) == 0)
                        vecindario.push_back(std::pair<size_t, size_t>(i, l));
                    C[i] = aux;
                }
            }
        }

        // Exploración aleatoria del entorno
        shuffle(vecindario.begin(), vecindario.end(), std::default_random_engine(seed+ev));

        // Exploración del entorno
        for (auto v = vecindario.begin(); v != vecindario.end() && !hay_mejora; ++v) {
            S = C;
            S[v->first] = v->second;
            f_vecino = fitness(S, X, R, clusters, lambda);
            if (f_vecino < f_actual) { // selección del primer mejor vecino
                f_actual = f_vecino;
                hay_mejora = true;
                C = S;
            }
        }
    }
    return C;
}

// ------------------------------ALGORITMOS GENÉTICOS------------------------------

// Variables generales
const size_t cromosomas = 50;
const float pc_agg = 0.7f, num_pm = 0.1f; // Probabilidad de cruce para AGG y numerador de probabilidad de mutación

// Funciones generales
int_matrix inicializar_poblacion(size_t n, size_t k) {
    int_matrix poblacion;
    std::vector<int> C;
    for (size_t p = 0; p < cromosomas; ++p) {
        while (empty_clusters(C, k) > 0) {
            C.clear();
            for (size_t i = 0; i < n; ++i)
                C.push_back(rand() % k);
        }
        poblacion.push_back(C);
    }
    return poblacion;
}

std::vector<double> evaluar_poblacion(const int_matrix& poblacion, const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda) {
    std::vector<double> evaluacion;
    for (auto& s : poblacion)
        evaluacion.push_back(fitness(s, X, R, clusters, lambda));
    return evaluacion;
}

int_matrix seleccion_generacional(const int_matrix& poblacion, const std::vector<double>& evaluacion) {
    int_matrix padres;
    size_t index1, index2;
    while (padres.size() < poblacion.size()) {
        index1 = rand()%poblacion.size();
        index2 = rand()%poblacion.size();
        padres.push_back(poblacion[evaluacion[index2] >= evaluacion[index1] ? index1 : index2]);
    }
    return padres;
}

std::vector<int> hijo_uniforme(const std::vector<int>& padre1, const std::vector<int>& padre2, std::vector<size_t> indices, size_t seed) {
    std::vector<int> hijo;
    size_t n = padre1.size();
    hijo.resize(n);
    shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));
    while (!indices.empty()) {
        if (indices.size() > n*0.5)
            hijo[indices.back()] = padre1[indices.back()];
        else
            hijo[indices.back()] = padre2[indices.back()];
        indices.pop_back();
    }
    return hijo;
}

void cruce_uniforme(int_matrix& padres, size_t seed) {
    size_t n_cruces = pc_agg*padres.size()*0.5f, n = padres[0].size();
    std::vector<size_t> indices;
    for (size_t i = 0; i < n; ++i)
        indices.push_back(i);
    for (size_t i = 1; i < n_cruces; ++i) {
        // primer hijo
        padres[i-1] = hijo_uniforme(padres[i-1], padres[i], indices, seed+i);
        // segundo hijo
        padres[i] = hijo_uniforme(padres[i-1], padres[i], indices, seed+i);
    }
}

void mutacion_uniforme(int_matrix& intermedia, size_t k) {
    size_t n = intermedia[0].size();
    size_t n_mutaciones = num_pm/n*intermedia.size();
    for (size_t i = 0; i < n_mutaciones; ++i)
        intermedia[i][rand()%n] = rand()%k;
}

// Algoritmos
std::vector<int> AGG_UN(const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t seed)
{
    std::vector<int> mejor;
    int_matrix poblacion = inicializar_poblacion(X.size(), clusters.size()), seleccionados;
    std::vector<double> evaluacion = evaluar_poblacion(poblacion, X, R, clusters, lambda);
    for (size_t ev = 0; ev < 100000; ev += cromosomas) {
        seleccionados = seleccion_generacional(poblacion, evaluacion);
        cruce_uniforme(seleccionados, seed);
        mutacion_uniforme(seleccionados, clusters.size());
    }
    return mejor;
}

#endif //PRACTICA1_MH_PAR_ALGORITMOS_H