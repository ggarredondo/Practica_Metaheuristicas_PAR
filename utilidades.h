#ifndef PRACTICA1_MH_PAR_UTILIDADES_H
#define PRACTICA1_MH_PAR_UTILIDADES_H
#include <vector>
#include <list>
#include <fstream>
#include <functional>
#include <algorithm>

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
    return datos;
}

typedef std::vector<std::list<int> > R_matrix;

R_matrix matriz_a_lista(const int_matrix& m) {
    R_matrix result;
    for (auto& v : m)
        result.push_back(std::list<int>(v.begin(), v.end()));
    return result;
}

class cluster {
private:
    std::vector<double> centroide;
    double_matrix puntos;

public:
    cluster(std::vector<double>& c) : centroide(c) {}

    inline void aniadir_punto(const std::vector<double>& punto) {
        puntos.push_back(punto);
    }

    void actualizar_centroide() {
        centroide = puntos[0];
        size_t size = puntos.size();
        for (size_t i = 1; i < size; ++i)
            transform(centroide.begin(), centroide.end(), puntos[i].begin(), centroide.begin(), std::plus<double>());
        for (size_t i = 0; i < size; ++i)
            centroide[i] /= size;
    }

    inline bool operator==(const cluster& c) const {
        return centroide == c.centroide;
    }
};

#endif //PRACTICA1_MH_PAR_UTILIDADES_H