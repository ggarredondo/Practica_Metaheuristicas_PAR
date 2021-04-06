#ifndef PRACTICA1_MH_PAR_ALGORITMOS_H
#define PRACTICA1_MH_ALGORITMOS_H
#include <random>
#include "utilidades.h"

void greedy_copkm(const double_matrix& X, const int_matrix& R, const R_list& Rl, std::vector<cluster>& C, size_t seed)
{
    std::vector<size_t> vRSI;
    for (size_t i = 0; i < X.size(); ++i)
        vRSI.push_back(i);
    shuffle(vRSI.begin(), vRSI.end(), std::default_random_engine(seed));
    std::list<size_t> RSI(vRSI.begin(), vRSI.end());

    do {
        for (auto& i : RSI) {
            
        }
    } while(true);
    
}

#endif //PRACTICA1_MH_PAR_ALGORITMOS_H