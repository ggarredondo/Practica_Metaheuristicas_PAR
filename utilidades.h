#ifndef PRACTICA1_MH_PAR_UTILIDADES_H
#define PRACTICA1_MH_PAR_UTILIDADES_H
#include <vector>
#include <list>
#include <fstream>
#include <functional>
#include <algorithm>
#include <cmath>

typedef std::vector<std::vector<int> > int_matrix;
typedef std::vector<std::vector<double> > double_matrix;

template <class T>
std::vector<std::vector<T> > archivo_a_matriz(std::ifstream archivo)
{
    std::vector<std::vector<T> > datos;
    std::vector<T> aux;
    std::string linea, number;
    while (archivo >> linea) {
        for (auto& i : linea) {
            if (i == ',') {
                aux.push_back(stod(number));
                number.clear();
            }
            else
                number.push_back(i);
        }
        aux.push_back(stod(number));
        number.clear();
        datos.push_back(aux);
        aux.clear();
    }
    archivo.close();
    return datos;
}

typedef std::vector<std::list<std::pair<size_t, int> > > R_matrix;

R_matrix matriz_a_lista(const int_matrix& m) {
    R_matrix result;
    std::list<std::pair<size_t, int> > aux;
    for (size_t i = 0; i < m.size(); ++i) {
        for (size_t j = 0; j < m[0].size(); ++j) {
            if (m[i][j] != 0 && i != j)
                aux.push_back(std::pair<size_t, int>(j, m[i][j]));
        }
        result.push_back(aux);
        aux.clear();
    }
    return result;
}

inline double distancia_euclidea(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0;
    for (size_t i = 0; i < a.size(); ++i)
        result += (a[i]-b[i])*(a[i]-b[i]);
    return sqrt(result);
}

inline double gen() {
    return (rand()%1000)*0.001;
}

// Pre: todos los puntos deben ser del mismo número de elementos
class cluster {
private:
    std::vector<double> centroide;
    double_matrix puntos;

public:
    cluster() {}

    cluster(size_t seed, size_t n) {
        centroide.resize(n);
        generate(centroide.begin(), centroide.end(), gen);
    }

    inline void aniadir_punto(const std::vector<double>& punto) {
        puntos.push_back(punto);
    }

    inline void eliminar_puntos() {
        puntos.clear();
    }

    void actualizar_centroide() {
        if (!puntos.empty()) {
            centroide = puntos[0];
            for (size_t i = 1; i < puntos.size(); ++i)
                transform(centroide.begin(), centroide.end(), puntos[i].begin(), centroide.begin(), std::plus<double>());
            size_t size = centroide.size();
            for (size_t i = 0; i < size; ++i)
                centroide[i] /= size;
        }
    }

    inline double distancia_centroide(const std::vector<double>& punto) const {
        return distancia_euclidea(centroide, punto);
    }

    void actualizar_puntos(const std::vector<int>& C, int ci, const double_matrix& X) {
        puntos.clear();
        for (size_t i = 0; i < C.size(); ++i) {
            if (ci == C[i])
                puntos.push_back(X[i]);
        }
        actualizar_centroide();
    }

    double distancia_intracluster(const std::vector<int>& C, int ci, const double_matrix& X) {
        actualizar_puntos(C, ci, X);
        double result = 0;
        for (auto &xi : puntos)
            result += distancia_centroide(xi);
        return result/puntos.size();
    }
};

inline size_t empty_clusters(const std::vector<int>& C, size_t k) {
    size_t sum = 0;
    for (int i = 0; i < k; ++i)
        sum += (count(C.begin(), C.end(), i)==0);
    return sum;
}

double desviacion_general(const std::vector<int>& C, const double_matrix& X, std::vector<cluster>& clusters) {
    double result = 0;
    size_t size = clusters.size();
    for (size_t i = 0; i < size; ++i)
        result += clusters[i].distancia_intracluster(C, i, X);
    return result/size;
}

double distancia_maxima(const double_matrix& X) {
    double max = 0, d;
    for (size_t i = 0; i < X.size(); ++i) {
        for (size_t j = i+1; j < X[0].size(); ++j) {
            d = distancia_euclidea(X[i], X[j]);
            if (d > max)
                max = d;
        }
    }
    return max;
}

#endif //PRACTICA1_MH_PAR_UTILIDADES_H