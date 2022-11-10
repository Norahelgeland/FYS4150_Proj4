

#include "Ising_model.hpp"
#include <armadillo>
#include <iostream>
#include <math.h>

std::mt19937 generator (123);
std::uniform_int_distribution<int> disINT(0.0, 19.);
std::uniform_real_distribution<double> disDOUBLE(0.0, 1.0);


Ising_model::Ising_model(double T_in, int L_in){
    

    //double random = dis(generator);

    // Setting values
    T = T_in;
    L = L_in;
    N = L*L;
    //S = arma::mat(L+2, L+2).fill(1);
    S = arma::mat(L, L).fill(1);
    Z = Z_fun();
    
}
    // HOW (???)
    double Ising_model::tot_energy(arma::mat S){
        double E_tot = 0;

        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                E_tot += S(i,j)*(S((i+1)%L,j));
                E_tot += S(i,j)*(S(i,(j+1)%L));

            }

         /* for(int i = 0; i < L-1; i++){
            for(int entry = -1; entry <= L; entry++){
                //E_tot += S.row(i)(entry) * S.row(i+1)(entry);
                //E_tot += S.col(i)(entry) * S.col(i+1)(entry);

                E_tot += S.row(i)(entry%L) * S.row(i+1)(entry%L);
                E_tot += S.col(i)(entry%L) * S.col(i+1)(entry%L);
            }  */
         }
        

/*          for(int i = 1; i <= L; i++){
            for(int j = 1; j <= L; j++){
                E_tot += S(i,j)*(S(i-1,j)+S(i,j-1)+S(i+1,j)+S(i,j+1));
            }
        }  */ 
        return -E_tot;
    }



    arma::vec Ising_model::possible_E(){
        arma::mat A = arma::mat(L, L).fill(1);

        arma::vec E_poss = arma::vec(N+1).fill(0);

        E_poss(0) = tot_energy(A);
        
        int counter = 1;
        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                A(i,j) = -1;
                E_poss(counter) = tot_energy(A);
                counter += 1;
            }
        } 

        return E_poss;
    }

    double Ising_model::tot_magnetization(arma::mat S){
        int Msum = 0;
        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                Msum += S(i,j);
            } 
        } 
    return Msum;
    }


    arma::vec Ising_model::possible_M(){

        arma::mat A = arma::mat(L, L).fill(1);

        arma::vec M_poss = arma::vec(N+1).fill(0);

        M_poss(0) = tot_magnetization(A);
        
        int counter = 1;
        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                A(i,j) = -1;
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
        double distr = (1./Z)*exp(-(1/T)*tot_energy(S)); 
        return distr;
    }

/*     arma::vec Ising_model::Possible_p(){
        
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
    } */

/*     // Correct with times pdf or mean (???)
    double Ising_model::Exp_value(arma::vec input){
        arma::vec pdf = Possible_p();
        double exp_val = 0;

        for(int s = 0; s < N+1; s++){
            exp_val += input(s)*pdf(s);
        }

        return exp_val;
    } */

    void Ising_model::MCMC(){

        arma::mat S_i = S;
        int num_i = disINT(generator);
        int num_j = disINT(generator);
         
        // Pick random spinn position and turn the spinn
        int old_sum = (S((num_i-1)%L, num_j) + S(num_i, (num_j-1)%L) + S((num_i+1)%L, num_j) + S(num_i, (num_j+1)%L));
        S(num_i, num_j) = S(num_i, num_j)*(-1);
        int sum_neighbours = (S((num_i-1)%L, num_j) + S(num_i, (num_j-1)%L) + S((num_i+1)%L, num_j) + S(num_i, (num_j+1)%L));
        
        double delta_E = sum_neighbours - old_sum;

        double r = disDOUBLE(generator);
        double p = exp(-(1/T)*delta_E);
        //double p = boltzmann_dist(S)/boltzmann_dist(S_i);
        //std::cout << p << " ";

        if(p>1){
            p = 1;
        }

        if(r > p){
            S = S_i;
        }

    }

    void Ising_model::update(){

            double E_col = 0;
            double B_col = 0;

            for(int j = 0; j < N; j++){
                MCMC();
                E_col += tot_energy(S);
                B_col += tot_magnetization(S);
                }  
            //std::cout << std::endl;

            // Update expectation values
            exp_val_E = E_col/N;
            exp_val_e = E_col/(N*N);
            exp_val_m = abs(B_col)/(N*N);
            spes_heat = 1./N*1./pow(T,2)*(((E_col*E_col)/N)-(exp_val_E*exp_val_E));
            suscept = 1./N*1./T*((abs(B_col*B_col)/(N*N))-(exp_val_m*exp_val_m));

    }
    






    



    


