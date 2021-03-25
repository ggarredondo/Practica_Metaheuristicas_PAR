#ifndef PRACTICA1_MH_PAR_PUNTO_H
#define PRACTICA1_MH_PAR_PUNTO_H
#include <vector>
#include <cmath>

class Punto {
private:
    std::vector<double> atributos;
    size_t cluster; // 0 == ningún cluster, 1 == cluster 1, 2 == cluster 2...

public:
    Punto(std::vector<double>& v) : atributos(v), cluster(0) {}

    inline void setCluster(size_t cluster) {
        this->cluster = cluster;
    }

    double distancia(const Punto& p) {
        double result = 0;
        for (size_t i = 0; i < atributos.size(); ++i)
            result += (this->atributos[i] - p.atributos[i])*(this->atributos[i] - p.atributos[i]);
        return sqrt(result);
    }
};


#endif //PRACTICA1_MH_PAR_PUNTO_H
