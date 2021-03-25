#include <iostream>
#include <fstream>
#include "punto.h"

std::vector<std::vector<int> > leer_const(std::ifstream archivo)
{
    std::vector<std::vector<int> > result;
    std::vector<int> aux;
    std::string linea;
    bool test = true;
    while (archivo >> linea) {
        if (test) {
            std::cout << linea << std::endl;
            test = false;
        }
    }
    return result;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    leer_const(std::ifstream("data/bupa_set_const_10.const"));
    return 0;
}
