#include <iostream>
#include <functional>
#include <vector>
#include <array>

#include "rscg.cpp"

using namespace std;

template <int N>
vector <double> matmul(const std::array<std::array <double, N>,N> &H, vector <double> &x){
    vector <double> Ax(N);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
        Ax[i] += H[i][j]*x[j];
        };
    };
    return Ax;
};

int main() {
    std::array<std::array<double, 3>, 3> H{ {
        { 1, 0, 3 },
        { 0, 2, 0 },
        { 3, 0, 5 }
    } };
    vector <double> x(3,1.0);
    auto f = std::bind(matmul<3>,H,std::placeholders::_1);

    RSCG* rscg;
    rscg = new RSCG(3,f);
    vector<complex <double> > vec_z(2,0.2);
    vec_z[0] = 0.2*1i;
    vec_z[1] = 4.0 + 0.5*1i;
    cout << vec_z[0] << std::endl;
    vector<complex <double> > Gij = rscg->calc_greenfunction_realH(0,0,vec_z);
    cout << "RSCG" << Gij[0] << " exact: 1.12377 - 0.383299 I" << std::endl;
    cout << "RSCG" << Gij[1] << " exact: 0.0844022 - 0.0339264 I" <<std::endl;

    std::array<std::array<double, 4>, 4> H2{
        {
            {0.998383,  1.52116,   1.49684,   0.556466},
            {1.52116,   0.425983,  0.164144,  0.571032},
            {1.49684,   0.164144,  0.512426,  0.806449},
            {0.556466,  0.571032,  0.806449,  1.22255}
        }
    };
    
    auto f2 = std::bind(matmul<4>,H2,std::placeholders::_1);
    rscg -> set_matvec(4,f2);
    vec_z[0] = 2.0 - 3.0*1i;
    vec_z[1] = -2.1 + 0.9*1i;
    cout << vec_z[0] << std::endl;
    Gij = rscg->calc_greenfunction_realH(0,0,vec_z);
    cout << "RSCG" << Gij[0] << " exact: 0.03436+0.221355 I" << std::endl;
    cout << "RSCG" << Gij[1] << " exact: -0.35166-0.307711 I" <<std::endl;

};