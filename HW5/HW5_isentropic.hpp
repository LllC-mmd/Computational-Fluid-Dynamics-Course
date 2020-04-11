#include <algorithm>

template<typename _Tv>
class Matrix{        
    public:
        _Tv *_elems;
        int nx;
        int ny;

        Matrix(int const& nj, int const& ni){
            ny = nj;
            nx = ni;
            _elems = new _Tv[nj * ni];
        }

        Matrix(const Matrix<_Tv>& sourceMat){
            ny = sourceMat.ny;
            nx = sourceMat.nx;
            _elems = new _Tv[ny * nx];
            for (int i = 0; i < ny; i++){
                for (int j = 0; j < nx; j++){
                    _elems[i*nx+j] = sourceMat[i][j];
                }
            }
        }

        friend void swap(Matrix<_Tv>& first, Matrix<_Tv>& second){
            using std::swap;
            swap(first.nx, second.nx);
            swap(first.ny, second.ny);
            swap(first._elems, second._elems);
        }

        Matrix& operator=(Matrix<_Tv> rhs){
            swap(*this, rhs);
            return *this;
        }

        _Tv* const* operator[](int const& j) const{
            return &_elems[j*nx];
        }

        _Tv* operator[](int const& j){
            return &_elems[j*nx];
        }

        ~Matrix(){
            delete[] _elems;
        }
};

struct Var{
    double v1;
    double v2;
    double v3;
    double v4;
};

