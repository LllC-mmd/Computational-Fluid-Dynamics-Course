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