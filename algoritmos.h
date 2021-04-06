#ifndef PRACTICA1_MH_PAR_ALGORITMOS_H
#define PRACTICA1_MH_PAR_ALGORITMOS_H
#include <random>
#include "utilidades.h"

// cj es igual al cluster al que pertenece xk y la restricción es cannot-link
// o
// cj no es igual al cluster al que pertenece xk y la restrición es must-link
size_t infeasibility(size_t xi, size_t cj, const std::vector<size_t> C, const R_matrix& R) {
    size_t xk = 0, inf = 0;
    for (auto &v : R[xi]) {
        if (xi != xk)
            inf += (cj == C[xk] && v == -1) || (cj != C[xk] && v == 1);
        xk++;
    }
    return inf;
}

// Para C, 0 == ningún cluster, 1 == cluster 1, 2 == cluster 2....
std::vector<size_t> greedy_copkm(const double_matrix& X, const R_matrix& R, std::vector<cluster>& clusters, size_t seed)
{
    std::vector<size_t> vRSI, C(X.size());
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