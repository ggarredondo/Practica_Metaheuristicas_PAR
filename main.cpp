#include <iostream>
#include "algoritmos.h"

int main() {

    int_matrix m = { {1,-1,1,1},
                     {-1,1,0,0},
                     {1,0,1,0},
                     {1,0,0,1} };

    for (size_t i = 0; i < m.size(); ++i) {
        for (size_t j = 0; j < m[i].size(); j++)
            std::cout << m[i][j] <<" ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    R_matrix R = matriz_a_lista(m);

    for (size_t i = 0; i < m.size(); ++i) {
        for (auto &j : m[i])
            std::cout << j <<" ";
        std::cout << std::endl;
    }

    size_t sum = 0;
    sum += false;
    std::cout << sum << std::endl;
}
