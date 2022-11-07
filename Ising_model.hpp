#ifndef __filename_hpp__
#define __filename_hpp__
#include <armadillo>

class Ising_model{

    private:

    public:
        double T;
        int L;
        int N;
        arma::mat S;
        int Z;

        Ising_model(){};
        // Constructor
        Ising_model(double T, int L);


        double tot_energy(arma::mat S);
        
        arma::vec possible_E();

        double Z_fun();

        double tot_magnetization(arma::mat S);

        arma::vec possible_M();

        double boltzmann_dist(arma::mat S);

        arma::vec Possible_p();
        
        double Exp_value(arma::vec input);

        void MCMC();

};

#endif
