#include <iostream>
#include "utilidades.h"

int main() {
    int_matrix matriz_r = archivo_a_matriz<int>(std::ifstream("data/bupa_set_const_10.const"));
    R_list lista_r = matriz_a_lista(matriz_r);
    size_t counter = 0;
    for (auto i = lista_r.begin(); i != lista_r.end() && counter < 344; ++i) {
        counter++;
        std::cout << i->r << ", ";
    }
    return 0;
}
