

#include "Ising_model.hpp"
#include <armadillo>
#include <iostream>
#include <math.h>

Ising_model::Ising_model(double T_in, int L_in){

    // Setting values
    T = T_in;
    L = L_in;
    N = L*L;
    S = arma::mat(L+2, L+2).fill(1);
    Z = Z_fun();
    
  /*   exp_val_E = Exp_value(possible_E());
    exp_val_e = Exp_value(possible_E()/N);
    exp_val_m = Exp_value(abs(possible_M()/N));
    spes_heat = 1./N*1./pow(T,2)*(Exp_value(pow(possible_E(),2))-pow(Exp_value(possible_E()),2));
    suscept = 1./N*1./T*(Exp_value(pow(possible_M(),2))-pow(Exp_value(abs(possible_M())),2)); */

}

    double Ising_model::tot_energy(arma::mat S){
        double E_tot = 0;

         for(int i = 0; i <= L; i++){
            for(int entry = 1; entry <= L; entry++){
            E_tot += S.row(i)(entry) * S.row(i+1)(entry);
            E_tot += S.col(i)(entry) * S.col(i+1)(entry);
            } 
         }
        

/*          for(int i = 1; i <= L; i++){
            for(int j = 1; j <= L; j++){
                E_tot += S(i,j)*(S(i-1,j)+S(i,j-1)+S(i+1,j)+S(i,j+1));
            }
        }  */ 
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

    double Ising_model::tot_magnetization(arma::mat S){
        int Msum = 0;
        for(int i = 1; i <= L; i++){
            for(int j = 1; j <= L; j++){
                Msum += S(i,j);
            } 
        } 
    return Msum;
    }


    arma::vec Ising_model::possible_M(){
        arma::mat A = arma::mat(L+2, L+2).fill(1);

        arma::vec M_poss = arma::vec(N+1).fill(0);

        M_poss(0) = tot_magnetization(A);
        
        int counter = 1;
        for(int i = 1; i <= L; i++){
            for(int j = 1; j <= L; j++){
                A(i,j) = -1;
                A.col(0) = A.col(L);
                A.col(L+1) = A.col(1);
                A.row(0) = A.row(L);
                A.row(L+1) = A.row(1);
                M_poss(counter) = tot_magnetization(A);
                counter += 1;
            }
        } 

        return M_poss;
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
    double Ising_model::Exp_value(arma::vec input){
        arma::vec pdf = Possible_p();
        double exp_val = 0;

        for(int s = 0; s < N+1; s++){
            exp_val += input(s)*pdf(s);
        }

        return exp_val;
    }

    void Ising_model::MCMC(){

        arma::mat S_i = S;
        int range = L - 1 + 1;
        int num_i = rand() % range + 1; 

        int num_j = rand() % range + 1; 

        S(num_i, num_j) = S(num_i, num_j)*(-1);
        S.col(0) = S.col(L);
        S.col(L+1) = S.col(1);
        S.row(0) = S.row(L);
        S.row(L+1) = S.row(1);
        double r = (rand() % 10)/10.;
        double p = boltzmann_dist(S)/boltzmann_dist(S_i);

        arma::vec E_col = arma::vec(N+1).fill(0);
        

        if(p>1){

            p = 1;
        }

        if(r > p){
            S = S_i;
        }

    }

    void Ising_model::update(){

            arma::vec E_col = arma::vec(N+1).fill(0);
            arma::vec B_col = arma::vec(N+1).fill(0);
            for(int j = 0; j < N; j++){
                MCMC();
                E_col(j) = tot_energy(S);
                B_col(j) = tot_magnetization(S);
                }   
         // Update expectation values

            exp_val_E = accu(E_col)/N;
            exp_val_e = accu(E_col/(N))/(N);
            exp_val_m = accu(abs(abs(B_col)/(N))/(N));
            //spes_heat
            //suscept

            /* exp_val_E = Exp_value(E_col);
            exp_val_e = Exp_value(E_col/N);
            exp_val_m = Exp_value(abs(B_col/N));
            spes_heat = 1./N*1./pow(T,2)*(Exp_value(pow(E_col,2))-pow(Exp_value(E_col),2));
            suscept = 1./N*1./T*(Exp_value(pow(B_col,2))-pow(Exp_value(abs(B_col)),2)); */
    }
    






    



    


