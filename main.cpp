#include <iostream>
#include "utilidades.h"

int main() {
    int_matrix test = archivo_a_matriz<int>(std::ifstream("data/bupa_set_const_10.const"));
    for (size_t i = 0; i < test.size(); ++i) {
        for (size_t j = 0; j < test[0].size(); ++j)
            std::cout << test[i][j] << " ";
        break;
        std::cout << std::endl;
    }
    return 0;
}
