#include <iostream>
#include "algoritmos.h"

int main() {
    size_t seed = 1;
    srand(seed);
    double_matrix data = archivo_a_matriz<double>(std::ifstream("data/bupa_set.dat"));

    std::vector<double> centroide = {1,1,1,1};
    std::vector<double> test = {2,3,4,5};
    std::vector<double> puntos = {2,3,4,3};
    transform(centroide.begin(), centroide.end(), puntos.begin(), centroide.begin(), std::plus<int>());
    std::list<size_t> lista = {1,2,3,4};
    for (auto& i : lista)
        std::cout << i << std::endl;

    if (test==puntos)
        std::cout << "yes" << std::endl;
    for (auto &v : centroide)
        std::cout << v << " ";
}
