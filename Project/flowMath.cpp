#include "flowMath.hpp"


double sgn(double x){
    if (x>=0){
        return 1;
    }
    else{
        return 0;
    }
}


double dReLU(double x){
    if (x>=0){
        return 1.0;
    }
    else{
        return 0.0;
    }
}

// ReFac(m, p) = m * (m-1) * ... * (m-p+1) 
int ReFac(int m, int p){
    int fac = 1;
    int c = std::min(m, p);

    for (int i = 0; i < c; i++){
        fac *= (m-i);
    }
    return fac;
}
