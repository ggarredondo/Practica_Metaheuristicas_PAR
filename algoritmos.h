#ifndef PRACTICA1_MH_PAR_ALGORITMOS_H
#define PRACTICA1_MH_PAR_ALGORITMOS_H
#include <random>
#include <cfloat>
#include "utilidades.h"
#include <chrono>
#include <cmath>

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
    std::vector<size_t> vRSI(X.size());
    iota(vRSI.begin(), vRSI.end(), 0);

    // Bajaramos los índices para recorrer X de forma aleatoria y sin repetición
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

// ------------------------------------------------------------------------------------

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

std::vector<int> generar_solucion_aleatoria(size_t n, size_t k) {
    std::vector<int> C(n);
    while (empty_clusters(C, k) > 0) {
        C.clear();
        for (size_t i = 0; i < n; ++i)
            C.push_back(rand()%k);
    }
    return C;
}

const size_t max_evaluaciones = 100000;

// ------------------------------BÚSQUEDAS BASADAS EN TRAYECTORIAS------------------------------

std::vector<int> busqueda_local(const std::vector<int>& ini, const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t max_ev, size_t seed)
{
    std::vector<std::pair<size_t, size_t> > vecindario;
    size_t k = clusters.size(), n = X.size(), aux;
    std::vector<int> C = ini, S;

    double f_actual = fitness(C, X, R, clusters, lambda), f_vecino;
    bool hay_mejora = true;
    for (size_t ev = 0; ev < max_ev && hay_mejora;) {
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
            ev++;
            if (f_vecino < f_actual) { // selección del primer mejor vecino
                f_actual = f_vecino;
                hay_mejora = true;
                C = S;
            }
        }
    }
    return C;
}

// Enfriamiento Simulado

inline double beta(double T0, double Tf, double M) {
    return (T0-Tf)/(M*T0*Tf);
}

inline double esquema_cauchy(double Tk, double T0, double Tf, double M) {
    return Tk/(1 + beta(T0,Tf,M)*Tk);
}

std::vector<int> cambio_cluster(const std::vector<int>& S, size_t n, size_t k) {
    std::vector<int> Si;
    do {
        Si = S;
        Si[rand() % n] = rand() % k;
    } while (S == Si || empty_clusters(Si, k) > 0);
    return Si;
}

std::vector<int> enfriamiento_simulado(const std::vector<int>& ini, const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t max_ev, size_t seed)
{
    size_t n = X.size(), k = clusters.size(), max_vecinos = 10*n, exitos = 1, max_exitos = 0.1*max_vecinos, M = max_ev/max_vecinos;
    std::vector<int> S = ini, Si, mejor_S = S;
    double fs = fitness(S, X, R, clusters, lambda), T0 = 0.3*fs/(-log(0.3)), Tk = T0, Tf = 0.003, fi, mejor_f = fs, dif_f;
    for (size_t ev = 0; ev < max_ev && exitos > 0;) {
        exitos = 0;
        for (size_t vecinos = 0; vecinos < max_vecinos && exitos < max_exitos; ++vecinos) {
            Si = cambio_cluster(S, n, k);
            fi = fitness(Si, X, R, clusters, lambda);
            ev++;
            dif_f = fi - fs;
            if (dif_f < 0 || rand()%10*0.1 < exp(-dif_f/Tk)) {
                exitos++;
                S = Si;
                fs = fi;
                if (fs < mejor_f) {
                    mejor_S = S;
                    mejor_f = fs;
                }
            }
        }
        Tk *= 0.9;
    }
    return mejor_S;
}

// Búsqueda multiarranque básica
std::vector<int> BMB(const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t seed)
{
    std::vector<int> S, mejor_S;
    size_t iteraciones = 10;
    double f, mejor_f = DBL_MAX;
    for (size_t ev = 0; ev < max_evaluaciones; ev += max_evaluaciones/iteraciones) {
        S = busqueda_local(generar_solucion_aleatoria(X.size(), clusters.size()), X, R, clusters, lambda, max_evaluaciones/iteraciones, seed+ev);
        f = fitness(S, X, R, clusters, lambda);
        if (f < mejor_f) {
            mejor_S = S;
            mejor_f = f;
        }
    }
    return mejor_S;
}

// Búsqueda Local Reiterada

void mutacion_ils(std::vector<int>& S, size_t tam, size_t k) {
    for (size_t i = rand()%S.size(), c = 0; c < tam; i = (i+1)%S.size(), ++c)
        S[i] = rand()%k;
}

std::vector<int> ILS(const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t seed)
{
    size_t iteraciones = 10;
    std::vector<int> S = busqueda_local(generar_solucion_aleatoria(X.size(), clusters.size()), X, R, clusters, lambda, max_evaluaciones/iteraciones, seed);
    std::vector<int> mejor_S = S;
    double f, mejor_f = DBL_MAX;
    for (size_t ev = max_evaluaciones/iteraciones; ev < max_evaluaciones; ev += max_evaluaciones/iteraciones) {
        mutacion_ils(S, 0.1*S.size(), clusters.size());
        S = busqueda_local(S, X, R, clusters, lambda, max_evaluaciones/iteraciones, seed+ev);
        f = fitness(S, X, R, clusters, lambda);
        if (f < mejor_f) {
            mejor_S = S;
            mejor_f = f;
        }
    }
    return mejor_S;
}

std::vector<int> ILS_ES(const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t seed)
{
    size_t iteraciones = 10;
    std::vector<int> S = enfriamiento_simulado(generar_solucion_aleatoria(X.size(), clusters.size()), X, R, clusters, lambda, max_evaluaciones/iteraciones, seed);
    std::vector<int> mejor_S = S;
    double f, mejor_f = DBL_MAX;
    for (size_t ev = max_evaluaciones/iteraciones; ev < max_evaluaciones; ev += max_evaluaciones/iteraciones) {
        mutacion_ils(S, 0.1*S.size(), clusters.size());
        S = enfriamiento_simulado(S, X, R, clusters, lambda, max_evaluaciones/iteraciones, seed+ev);
        f = fitness(S, X, R, clusters, lambda);
        if (f < mejor_f) {
            mejor_S = S;
            mejor_f = f;
        }
    }
    return mejor_S;
}

// ------------------------------ALGORITMOS GENÉTICOS------------------------------

//---Variables generales---
const size_t cromosomas = 50;
const float pc_agg = 0.7f, num_pm = 0.1f; // Probabilidad de cruce para AGG y numerador de probabilidad de mutación

//---Funciones generales---
int_matrix inicializar_poblacion(size_t n, size_t k) {
    int_matrix poblacion;
    for (size_t p = 0; p < cromosomas; ++p)
        poblacion.push_back(generar_solucion_aleatoria(n, k));
    return poblacion;
}

std::vector<double> evaluar_poblacion(const int_matrix& poblacion, const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda) {
    std::vector<double> evaluacion;
    for (auto& s : poblacion)
        evaluacion.push_back(fitness(s, X, R, clusters, lambda));
    return evaluacion;
}

//---Operadores de selección---
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

int_matrix seleccion_estacionario(const int_matrix& poblacion, const std::vector<double>& evaluacion) {
    int_matrix padres;
    size_t index1, index2;
    for (size_t i = 0; i < 2; ++i) {
        index1 = rand()%poblacion.size();
        index2 = rand()%poblacion.size();
        padres.push_back(poblacion[evaluacion[index2] >= evaluacion[index1] ? index1 : index2]);
    }
    return padres;
}

//---Operadores de cruce---

// Cruce uniforme
std::vector<int> hijo_uniforme(const std::vector<int>& padre1, const std::vector<int>& padre2, std::vector<size_t> indices, size_t seed) {
    std::vector<int> hijo;
    size_t n = padre1.size();
    hijo.resize(n);
    shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));
    while (!indices.empty()) {
        if (indices.size() > n*0.5f)
            hijo[indices.back()] = padre1[indices.back()];
        else
            hijo[indices.back()] = padre2[indices.back()];
        indices.pop_back();
    }
    return hijo;
}

void cruce_uniforme(float p, int_matrix& padres, size_t seed) {
    size_t n_cruces = p*padres.size()*0.5f;
    std::vector<size_t> indices(padres[0].size());
    iota(indices.begin(), indices.end(), 0);
    std::vector<int> aux;
    for (size_t i = 1; i < n_cruces*2; i += 2) {
        // primer hijo
        aux = padres[i-1];
        padres[i-1] = hijo_uniforme(padres[i-1], padres[i], indices, seed+i);
        // segundo hijo
        padres[i] = hijo_uniforme(aux, padres[i], indices, seed-i);
    }
}

// Cruce segmento fijo
std::vector<int> hijo_segmento_fijo(const std::vector<int>& padre1, const std::vector<int>& padre2, std::vector<size_t> indices, size_t seed) {
    std::vector<int> hijo;
    size_t n = padre1.size(), r = rand()%n, v = rand()%n, end = (r+v)%n;
    hijo.resize(n);
    for (size_t i = r; i != end; i = (i+1)%n) {
        hijo[i] = padre1[i];
        indices.erase(std::find(indices.begin(), indices.end(), i));
    }
    hijo[end] = padre1[end];
    indices.erase(std::find(indices.begin(), indices.end(), end));
    // resto de genes se seleccionan de manera uniforme
    n = indices.size();
    shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));
    while (!indices.empty()) {
        if (indices.size() > n*0.5f)
            hijo[indices.back()] = padre1[indices.back()];
        else
            hijo[indices.back()] = padre2[indices.back()];
        indices.pop_back();
    }
    return hijo;
}

void cruce_segmento_fijo(float p, int_matrix& padres, size_t seed) {
    size_t n_cruces = p*padres.size()*0.5f;
    std::vector<size_t> indices(padres[0].size());
    iota(indices.begin(), indices.end(), 0);
    std::vector<int> aux;
    for (size_t i = 1; i < n_cruces*2; i += 2) {
        // primer hijo
        aux = padres[i-1];
        padres[i-1] = hijo_segmento_fijo(padres[i-1], padres[i], indices, seed+i);
        // segundo hijo
        padres[i] = hijo_segmento_fijo(aux, padres[i], indices, seed-i);
    }
}

//---Operador de mutación---
void mutacion_uniforme(int_matrix& intermedia, size_t k) {
    size_t n = intermedia[0].size(), M = intermedia.size();
    float n_mutaciones = num_pm * M;
    for (size_t i = 0; i < n_mutaciones; ++i)
        intermedia[rand()%M][rand()%n] = rand()%k;
}

//---Algoritmos genéticos---
std::vector<int> AGG_UN(const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t seed)
{
    size_t index_mejor, index_peor;
    double ev_mejor;
    int_matrix poblacion = inicializar_poblacion(X.size(), clusters.size()), seleccionados;
    std::vector<double> evaluacion = evaluar_poblacion(poblacion, X, R, clusters, lambda);
    for (size_t ev = 0; ev < max_evaluaciones; ev += cromosomas) {
        // elitismo - encontrar el mejor de la población actual
        index_mejor = std::min_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin();
        ev_mejor = evaluacion[index_mejor];

        seleccionados = seleccion_generacional(poblacion, evaluacion);
        cruce_uniforme(pc_agg, seleccionados, seed+ev);
        mutacion_uniforme(seleccionados, clusters.size());
        evaluacion = evaluar_poblacion(seleccionados, X, R, clusters, lambda);

        // elitismo - reemplazamiento
        if (std::find(seleccionados.begin(), seleccionados.end(), poblacion[index_mejor]) == seleccionados.end()) {
           index_peor = std::max_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin();
           seleccionados[index_peor] = poblacion[index_mejor];
           evaluacion[index_peor] = ev_mejor;
        }
        poblacion = seleccionados;
    }
    return poblacion[std::min_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin()];
}

std::vector<int> AGG_SF(const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t seed)
{
    size_t index_mejor, index_peor;
    double ev_mejor;
    int_matrix poblacion = inicializar_poblacion(X.size(), clusters.size()), seleccionados;
    std::vector<double> evaluacion = evaluar_poblacion(poblacion, X, R, clusters, lambda);
    for (size_t ev = 0; ev < max_evaluaciones; ev += cromosomas) {
        // elitismo - encontrar el mejor de la población actual
        index_mejor = std::min_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin();
        ev_mejor = evaluacion[index_mejor];

        seleccionados = seleccion_generacional(poblacion, evaluacion);
        cruce_segmento_fijo(pc_agg, seleccionados, seed+ev);
        mutacion_uniforme(seleccionados, clusters.size());
        evaluacion = evaluar_poblacion(seleccionados, X, R, clusters, lambda);

        // elitismo - reemplazamiento
        if (std::find(seleccionados.begin(), seleccionados.end(), poblacion[index_mejor]) == seleccionados.end()) {
            index_peor = std::max_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin();
            seleccionados[index_peor] = poblacion[index_mejor];
            evaluacion[index_peor] = ev_mejor;
        }
        poblacion = seleccionados;
    }
    return poblacion[std::min_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin()];
}

std::vector<int> AGE_UN(const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t seed)
{
    size_t index_peor;
    int_matrix poblacion = inicializar_poblacion(X.size(), clusters.size()), seleccionados;
    std::vector<double> evaluacion = evaluar_poblacion(poblacion, X, R, clusters, lambda), ev_h;
    for (size_t ev = 0; ev < max_evaluaciones; ev += 2) {
        seleccionados = seleccion_estacionario(poblacion, evaluacion);
        cruce_uniforme(1, seleccionados, seed+ev);
        mutacion_uniforme(seleccionados, clusters.size());
        ev_h = evaluar_poblacion(seleccionados, X, R, clusters, lambda);

        // Reemplazamiento
        poblacion.push_back(seleccionados.front());
        evaluacion.push_back(ev_h.front());
        poblacion.push_back(seleccionados.back());
        evaluacion.push_back(ev_h.back());
        for (size_t i = 0; i < 2; ++i) {
            index_peor = std::max_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin();
            poblacion.erase(poblacion.begin()+index_peor);
            evaluacion.erase(evaluacion.begin()+index_peor);
        }
        ev_h.clear();
    }
    return poblacion[std::min_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin()];
}

std::vector<int> AGE_SF(const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t seed)
{
    size_t index_peor;
    int_matrix poblacion = inicializar_poblacion(X.size(), clusters.size()), seleccionados;
    std::vector<double> evaluacion = evaluar_poblacion(poblacion, X, R, clusters, lambda), ev_h;
    for (size_t ev = 0; ev < max_evaluaciones; ev += 2) {
        seleccionados = seleccion_estacionario(poblacion, evaluacion);
        cruce_segmento_fijo(1, seleccionados, seed+ev);
        mutacion_uniforme(seleccionados, clusters.size());
        ev_h = evaluar_poblacion(seleccionados, X, R, clusters, lambda);

        // Reemplazamiento
        poblacion.push_back(seleccionados.front());
        evaluacion.push_back(ev_h.front());
        poblacion.push_back(seleccionados.back());
        evaluacion.push_back(ev_h.back());
        for (size_t i = 0; i < 2; ++i) {
            index_peor = std::max_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin();
            poblacion.erase(poblacion.begin()+index_peor);
            evaluacion.erase(evaluacion.begin()+index_peor);
        }
        ev_h.clear();
    }
    return poblacion[std::min_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin()];
}

//---Algoritmo memético---

size_t busqueda_local_suave(std::vector<int>& cromosoma, double& evaluacion, const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, float epsilon, size_t seed) {
    size_t n_evaluaciones = 0, n = cromosoma.size(), fallos = 0, k = clusters.size();
    std::vector<int> S;
    std::vector<size_t> RSI(n);
    iota(RSI.begin(), RSI.end(), 0);
    shuffle(RSI.begin(), RSI.end(), std::default_random_engine(seed));
    double nueva_ev;
    bool mejora = true;
    for (size_t i = 0; i < n && (mejora || fallos < epsilon); ++i) {
        mejora = false;
        for (size_t l = 0; l < k; ++l) {
            if (cromosoma[RSI[i]] != l) {
                S = cromosoma;
                S[RSI[i]] = l;
                nueva_ev = fitness(S, X, R, clusters, lambda);
                ++n_evaluaciones;
                if (nueva_ev < evaluacion) {
                    mejora = true;
                    evaluacion = nueva_ev;
                    cromosoma = S;
                }
            }
        }
        fallos += !mejora;
    }
    return n_evaluaciones;
}

// pbl: probabilidad para aplicar BL
// mejor: si se aplica BLS a los pbl*N mejores de la población
std::vector<int> AM(const double_matrix& X, const R_list& R, std::vector<cluster>& clusters, double lambda, size_t seed, float pbl, bool mejor) {
    size_t index_mejor, index_peor;
    double ev_mejor;
    int_matrix poblacion = inicializar_poblacion(X.size(), clusters.size()), seleccionados;
    float n_explotaciones = pbl*poblacion.size();
    std::vector<double> evaluacion = evaluar_poblacion(poblacion, X, R, clusters, lambda);
    std::vector<size_t> indices(cromosomas);
    std::iota(indices.begin(), indices.end(), 0);

    for (size_t ev = 0, gen = 0; ev < max_evaluaciones; ev += cromosomas, ++gen) {
        // Cada 10 generaciones, aplicar BL
        if (gen%10 == 0 && gen > 0) {
            if (mejor)
                stable_sort(indices.begin(), indices.end(),[&evaluacion](size_t i1, size_t i2) {return evaluacion[i1] < evaluacion[i2];});
            for (size_t i = 0; i < n_explotaciones; ++i)
                ev += busqueda_local_suave(poblacion[indices[i]], evaluacion[indices[i]], X, R, clusters, lambda, 0.1f*X.size(), seed+ev+i);
        }

        // elitismo - encontrar el mejor de la población actual
        index_mejor = std::min_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin();
        ev_mejor = evaluacion[index_mejor];

        seleccionados = seleccion_generacional(poblacion, evaluacion);
        cruce_segmento_fijo(pc_agg, seleccionados, seed+ev);
        mutacion_uniforme(seleccionados, clusters.size());
        evaluacion = evaluar_poblacion(seleccionados, X, R, clusters, lambda);

        // elitismo - reemplazamiento
        if (std::find(seleccionados.begin(), seleccionados.end(), poblacion[index_mejor]) == seleccionados.end()) {
            index_peor = std::max_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin();
            seleccionados[index_peor] = poblacion[index_mejor];
            evaluacion[index_peor] = ev_mejor;
        }
        poblacion = seleccionados;
    }
    index_mejor = std::min_element(evaluacion.begin(), evaluacion.end()) - evaluacion.begin();
    return poblacion[index_mejor];
}

#endif //PRACTICA1_MH_PAR_ALGORITMOS_H