#ifndef PRACTICA1_MH_PAR_UTILIDADES_H
#define PRACTICA1_MH_PAR_UTILIDADES_H
#include <vector>
#include <fstream>

typedef std::vector<std::vector<int> > int_matrix;
typedef std::vector<std::vector<double> > double_matrix;

template <class T>
std::vector<std::vector<T> > archivo_a_matriz(std::ifstream archivo)
{
    std::vector<std::vector<T> > datos;
    std::vector<T> aux;
    std::string linea, number;
    while (archivo >> linea) {
        for (auto i = linea.begin(); i != linea.end(); ++i) {
            if (*i == ',') {
                aux.push_back(stod(number));
                number.clear();
            }
            else
                number.push_back(*i);
        }
        aux.push_back(stod(number));
        number.clear();
        datos.push_back(aux);
        aux.clear();
    }
    return datos;
}

#endif //PRACTICA1_MH_PAR_UTILIDADES_H