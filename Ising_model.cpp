

#include "Ising_model.hpp"
#include <armadillo>
#include <iostream>

Ising_model::Ising_model(double T_in, int L_in){

    // Setting values
    T = T_in;
    L = L_in;
    N = L*L;
    S = arma::mat(L+2, L+2).fill(1);
    Z = Z_fun();

}

    double Ising_model::tot_energy(arma::mat S){
        double E_tot = 0;
        for(int i = 1; i <= L; i++){
            for(int j = 1; j <= L; j++){
                E_tot += S(i,j)*(S(i-1,j)+S(i,j-1)+S(i+1,j)+S(i,j+1));
            }
        } 
        return -E_tot;
    }

    arma::vec Ising_model::possible_E(){
        arma::mat A = arma::mat(L+2, L+2).fill(1);

        arma::vec E_poss = arma::vec(N+1).fill(0);

        E_poss(0) = tot_energy(A);
        
        int counter = 1;
        for(int i = 1; i <= L; i++){
            for(int j = 1; j <= L; j++){
                A(i,j) = -1;
                A.col(0) = A.col(L);
                A.col(L+1) = A.col(1);
                A.row(0) = A.row(L);
                A.row(L+1) = A.row(1);
                E_poss(counter) = tot_energy(A);
                counter += 1;
            }
        } 

        return E_poss;
    }

    double Ising_model::tot_magnetization(){
        int Msum = 0;
        for(int i = 1; i <= L; i++){
            for(int j = 1; j <= L; j++){
                Msum += S(i,j);
            } 
        } 
    return Msum;
    }

    double Ising_model::Z_fun(){

        arma::vec Epos = possible_E();
        
         double Z = 0;
        for(int i = 0; i <= N; i++){
            Z += exp(-(1/T)*Epos(i));
        }
        return Z; 
    }

    double Ising_model::boltzmann_dist(arma::mat S){
        double distr = (1./Z)*exp(-(1/T)*tot_energy(S)); //(1/Z)*exp(-(1/T)*tot_energy(S));
        return distr;
    }

    arma::vec Ising_model::Possible_p(){
        
        arma::mat A = arma::mat(L+2, L+2).fill(1);
        arma::vec P_poss = arma::vec(N+1).fill(0);

        P_poss(0) = boltzmann_dist(A);
        
        int counter = 1;
        for(int i = 1; i <= L; i++){
            for(int j = 1; j <= L; j++){
                A(i,j) = -1;
                A.col(0) = A.col(L);
                A.col(L+1) = A.col(1);
                A.row(0) = A.row(L);
                A.row(L+1) = A.row(1);
                P_poss(counter) = boltzmann_dist(A);
                counter += 1;
            }
        } 

        return P_poss;
    }





    



    


