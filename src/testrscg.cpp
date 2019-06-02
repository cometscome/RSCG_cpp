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
    auto Gij = rscg->calc_greenfunction_realH(0,0,vec_z);
    cout << "RSCG" << Gij[0] << " exact: 1.12377 - 0.383299 I" << std::endl;
    cout << "RSCG" << Gij[1] << " exact: 0.0844022 - 0.0339264 I" <<std::endl;

    vector <int> vec_i(2,0);
    vec_i[1] = 2;
    auto mat_Gij = rscg->calc_greenfunction_realH<2>(vec_i,0,vec_z);
    cout << "RSCG" << mat_Gij[0][0] << " exact:   1.12377-0.383299 I" << std::endl;
    cout << "RSCG" << mat_Gij[1][0] << " exact:  -0.682371+0.202684 I" << std::endl;
    cout << "RSCG" << mat_Gij[0][1] << " exact:   0.0844022-0.0339264 I" << std::endl;
    cout << "RSCG" << mat_Gij[1][1] << " exact:  -0.243277-0.0198593 I" << std::endl; 
};