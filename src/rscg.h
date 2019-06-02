using namespace std;

#include <vector>
#include <cmath>
#include <memory>
#include <new>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <time.h>
#include <functional>
#include <array>

class RSCG{
    private:
    std::function<vector <complex <double> >(vector <complex <double> >&)> matvec;
    std::function<vector <double>(vector <double>&)> matvec_dble;
    int _n;
    int _maximumsteps; //Maximum iteration steps
    double _eps;
    

    public:
    RSCG(int n,std::function<vector <complex <double> >(vector <complex <double> >&)> f);
    RSCG(int n,std::function<vector <double>(vector <double>&)> f);

    vector <complex <double> > calc_greenfunction_complexH(int i,int j,
        vector <complex <double> > vec_z //z_j
        );

    template <int N> std::array<vector <complex <double> > ,N> calc_greenfunction_complexH(vector <int> vec_i,int j,
        vector <complex <double> > vec_z //z_j
        );

    vector <complex <double> > calc_greenfunction_realH(int i,int j,
        vector <complex <double> > vec_z //z_j
        );       

    template <int N> std::array<vector <complex <double> > ,N> calc_greenfunction_realH(vector <int> vec_i, int j,
        vector <complex <double> > vec_z //z_j
        );   


    void set_maximumsteps(int itemax){
        _maximumsteps = itemax;
    };

};