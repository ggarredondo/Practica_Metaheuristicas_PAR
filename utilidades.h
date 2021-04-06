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

struct R {
    size_t i, j;
    int r;
    R() : i(0), j(0), r(0) {}
    R(size_t i, size_t j, int r) : i(i), j(j), r(r) {}
};
typedef std::list<R> R_list;

R_list matriz_a_lista(int_matrix m) 
{
    R_list restricciones;
    R rest;
    for (size_t i = 0; i < m.size(); ++i) {
        for (size_t j = i+1; j < m[i].size(); ++j) {
            rest = R(i, j, m[i][j]);
            restricciones.push_back(rest);
        }
    }
    return restricciones;
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
};

#endif //PRACTICA1_MH_PAR_UTILIDADES_H