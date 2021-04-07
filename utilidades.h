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

typedef std::vector<std::list<int> > R_matrix;

R_matrix matriz_a_lista(const int_matrix& m) {
    R_matrix result;
    for (auto& v : m)
        result.push_back(std::list<int>(v.begin(), v.end()));
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
    size_t index;

public:
    cluster(size_t index, size_t seed, size_t size) {
        this->index = index;
        centroide.resize(size);
        generate(centroide.begin(), centroide.end(), gen);
    }

    inline void aniadir_punto(const std::vector<double>& punto) {
        puntos.push_back(punto);
    }

    void actualizar_centroide() {
        if (!puntos.empty()) {
            centroide = puntos[0];
            for (size_t i = 1; i < puntos.size(); ++i)
                transform(centroide.begin(), centroide.end(), puntos[i].begin(), centroide.begin(), std::plus<double>());
            size_t size = centroide.size();
            for (size_t i = 0; i < size; ++i)
                centroide[i] /= size;
            puntos.clear();
        }
    }

    inline double distancia_centroide(const std::vector<double>& punto) const {
        return distancia_euclidea(centroide, punto);
    }

    inline bool operator==(const cluster& c) const {
        return centroide == c.centroide;
    }

    inline bool empty(const std::vector<int>& C) const {
        return find(C.begin(), C.end(), index)==C.end();
    }
};

inline size_t empty_clusters(const std::vector<cluster>& clusters, const std::vector<int>& C) {
    size_t sum = 0;
    for (auto& ci : clusters)
        sum += ci.empty(C);
    return sum;
}

#endif //PRACTICA1_MH_PAR_UTILIDADES_H