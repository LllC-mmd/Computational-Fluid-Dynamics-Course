#include "variable.hpp"
#include "consts.hpp"


Var toVar(const U& u){
    Var var;

    var.v1 = u.h;
    var.v2 = u.hu;
    var.v3 = u.hv;

    return var;
}


Var toVar(const FluxF& f){
    Var var;

    var.v1 = f.F1;
    var.v2 = f.F2;
    var.v3 = f.F3;

    return var;
}


Var operator+(const Var& a, const Var& b){
    Var var_sum;

    var_sum.v1 = a.v1 + b.v1;
    var_sum.v2 = a.v2 + b.v2;
    var_sum.v3 = a.v3 + b.v3;

    return var_sum;
};


Var operator-(const Var& a, const Var& b){
    Var var_sub;

    var_sub.v1 = a.v1 - b.v1;
    var_sub.v2 = a.v2 - b.v2;
    var_sub.v3 = a.v3 - b.v3;

    return var_sub;
};


Var operator*(double alpha, const Var&a){
    Var var_mult;

    var_mult.v1 = alpha * a.v1;
    var_mult.v2 = alpha * a.v2;
    var_mult.v3 = alpha * a.v3;

    return var_mult;
};


Var operator/(const Var&a, double denominator){
    Var var_div;

    var_div.v1 = a.v1 / denominator;
    var_div.v2 = a.v2 / denominator;
    var_div.v3 = a.v3 / denominator;

    return var_div;
};


U operator+(const U& a, const U& b){
    U var_sum;

    var_sum.h = a.h + b.h;
    var_sum.hu = a.hu + b.hu;
    var_sum.hv = a.hv + b.hv;
    var_sum.eta = a.eta + b.eta;

    return var_sum;
}

U operator+(const U& a, const Var& b){
    U var_sum;

    var_sum.h = a.h + b.v1;
    var_sum.hu = a.hu + b.v2;
    var_sum.hv = a.hv + b.v3;
    var_sum.eta = a.eta;

    return var_sum;
}


U operator-(const U& a, const U& b){
    U var_sub;

    var_sub.h = a.h - b.h;
    var_sub.hu = a.hu - b.hu;
    var_sub.hv = a.hv - b.hv;
    var_sub.eta = a.eta - b.eta;

    return var_sub;
}

U operator*(double alpha, const U&a){
    U var_mult;

    var_mult.h = alpha * a.h;
    var_mult.hu = alpha * a.hu;
    var_mult.hv = alpha * a.hv;
    var_mult.eta = alpha * a.eta;

    return var_mult;
}

U operator/(const U&a, double denominator){
    U var_div;

    var_div.h = a.h / denominator;
    var_div.hu =  a.hu / denominator;
    var_div.hv = a.hv / denominator;
    var_div.eta = a.eta / denominator;

    return var_div;
}


FluxF operator+(const FluxF& a, const FluxF& b){
    FluxF var_sum;

    var_sum.F1 = a.F1 + b.F1;
    var_sum.F2 = a.F2 + b.F2;
    var_sum.F3 = a.F3 + b.F3;

    return var_sum;
}

FluxF operator-(const FluxF& a, const FluxF& b){
    FluxF var_sub;

    var_sub.F1 = a.F1 - b.F1;
    var_sub.F2 = a.F2 - b.F2;
    var_sub.F3 = a.F3 - b.F3;

    return var_sub;
}

FluxF operator*(double alpha, const FluxF& a){
    FluxF var_mult;

    var_mult.F1 = alpha * a.F1;
    var_mult.F2 = alpha * a.F2;
    var_mult.F3 = alpha * a.F3;

    return var_mult;
}

FluxF operator/(const FluxF& a, double denominator){
    FluxF var_div;

    var_div.F1 = a.F1 / denominator;
    var_div.F2 =  a.F2 / denominator;
    var_div.F3 = a.F3 / denominator;

    return var_div;
}

FluxF get_F(const U& u){
    FluxF f;

    double zb = u.eta - u.h;
    f.F1 = u.hu;
    f.F2 = ((u.h==0.0) ? 0.0 : u.hu*u.hu/u.h) + 0.5*g*(u.h*u.h-zb*zb);
    f.F3 = (u.h==0.0) ? 0.0 : u.hu*u.hv/u.h;

    return f;
}

FluxF toF(const Var& var){
    FluxF f;

    f.F1 = var.v1;
    f.F2 = var.v2;
    f.F3 = var.v3;

    return f;
}