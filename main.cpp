#include <iostream>
#include <fstream>
#include "punto.h"

std::vector<std::vector<int> > leer_const(std::ifstream archivo)
{
    std::vector<int> aux;
    std::string line;
    while (archivo >> line) {
        std::cout << line << std::endl;
    }
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    leer_const(std::ifstream("data/bupa_set_const_10.const"));
    return 0;
}
