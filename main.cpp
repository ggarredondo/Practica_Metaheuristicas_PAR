#include <iostream>
#include <fstream>
#include "punto.h"

typedef std::vector<std::vector<int> > int_matrix;

int_matrix leer_const(std::ifstream archivo)
{
    int_matrix restricciones;
    std::vector<int> aux;
    std::string linea, number;
    bool test = true;
    while (archivo >> linea) {
        for (auto i = linea.begin(); i != linea.end(); ++i) {
            if (*i == ',') {
                aux.push_back(stoi(number));
                number.clear();
            }
            else
                number.push_back(*i);
        }
        aux.push_back(stoi(number));
        restricciones.push_back(aux);
    }
    return restricciones;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    int_matrix test = leer_const(std::ifstream("data/bupa_set_const_10.const"));
    for (auto i = 0; i < test.size(); ++i) {
        for (auto j = 0; j < test.size(); ++j)
            std::cout << test[i][j] << " ";
        break;
    }
    return 0;
}
