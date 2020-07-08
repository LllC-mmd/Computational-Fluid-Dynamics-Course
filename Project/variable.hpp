#include <algorithm>

template<typename _Tv>
class Vector{        
    public:
        _Tv *_elems;
        int nx;

        Vector(int const& ni){
            nx = ni;
            _elems = new _Tv[ni];
        }

        Vector(const Vector<_Tv>& sourceVec){
            nx = sourceVec.nx;
            _elems = new _Tv[nx];
            for (int i = 0; i < nx; i++){
                _elems[i] = sourceVec[i];
            }
        }

        friend void swap(Vector<_Tv>& first, Vector<_Tv>& second){
            using std::swap;
            swap(first.nx, second.nx);
            swap(first._elems, second._elems);
        }

        Vector& operator=(Vector<_Tv> rhs){
            swap(*this, rhs);
            return *this;
        }

        const _Tv& operator[](int const& j) const{
            return _elems[j];
        }

        _Tv& operator[](int const& j){
            return _elems[j];
        }

        ~Vector(){
            delete[] _elems;
        }
};


struct Var{
    double v1;
    double v2;
    double v3;
};

// U[0] = h; U[1] = hu; U[2] = hv;
struct U{
    double h;
    double hu;
    double hv;
    double eta;
};

// FluxF[0] = h*u; FluxF[1] = h*u**2 + g*(h**2-zb**2)/2; FluxF[2] = h*u*v;
struct FluxF{
    double F1;
    double F2;
    double F3;
};

// S[0] = P; S[1] = -g*(h+zb)*i-g*h*Sfx; S[2] = -g*(h+zb)*i-g*h*Sfy
struct Source{
    double S1;
    double S2;
    double S3;
};

Var toVar(const U& u);
Var toVar(const FluxF& f);

Var operator+(const Var& a, const Var& b);
Var operator-(const Var& a, const Var& b);
Var operator*(double alpha, const Var&a);
Var operator/(const Var&a, double denominator);

U operator+(const U& a, const U& b);
U operator+(const U& a, const Var& b);
U operator-(const U& a, const U& b);
U operator*(double alpha, const U&a);
U operator/(const U&a, double denominator);

FluxF operator+(const FluxF& a, const FluxF& b);
FluxF operator-(const FluxF& a, const FluxF& b);
FluxF operator*(double alpha, const FluxF& a);
FluxF operator/(const FluxF& a, double denominator);
FluxF get_F(const U& u);
FluxF toF(const Var& var);