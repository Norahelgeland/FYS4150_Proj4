

#include "Ising_model.hpp"
#include <armadillo>
#include <iostream>
#include <math.h>

std::mt19937 generator (123);
std::uniform_int_distribution<int> disINT(0.0, 39.);
std::uniform_real_distribution<double> disDOUBLE(0.0, 1.0);


Ising_model::Ising_model(double T_in, int L_in){
    
    // Setting values
    T = T_in;
    L = L_in;
    N = L*L;
    S = arma::mat(L, L).fill(1);
    boltzmann_factor = make_boltzmann_factors();
    E_tot = tot_energy(S);
    
}
    double Ising_model::tot_energy(arma::mat S){
        double E_tot = 0;

        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                E_tot += S(i,j)*(S((i+1)%L,j));
                E_tot += S(i,j)*(S(i,(j+1)%L));
            }
         }
        // std::cout << S;
        // std::cout << "\n\n E tot:";
        // std::cout << -E_tot;
        // std::cout << "\n\n";
        
        return -E_tot;
    }

    std::map <double, double> Ising_model::make_boltzmann_factors(){
        std::map <double, double> boltzmann_factor;
        boltzmann_factor[-8.] = exp(-(1/T)*-8.);
        boltzmann_factor[-4.] = exp(-(1/T)*-4.);
        boltzmann_factor[0.] = exp(-(1/T)*0.);
        boltzmann_factor[4.] = exp(-(1/T)*4.);
        boltzmann_factor[8.] = exp(-(1/T)*8.);
        return boltzmann_factor;
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


/*     arma::vec Ising_model::possible_M(){

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
    } */

/* 
    double Ising_model::boltzmann_dist(arma::mat S){
        double distr = (1./Z)*exp(-(1/T)*tot_energy(S)); 
        return distr;
    } */

    void Ising_model::MCMC(){

        int num_i = disINT(generator);
        int num_j = disINT(generator);
/*         std::cout << "\n\n i:";
        std::cout << num_i;
        std::cout << "\n\n j:";
        std::cout << num_j; */
        // Pick random spinn position and turn the spinn
        int old_sum =  S(num_i, num_j)*(S((num_i-1+L)%L, num_j) + S(num_i, (num_j-1+L)%L) + S((num_i+1)%L, num_j) + S(num_i, (num_j+1)%L));
        //std::cout << "\n\n old_sum:";
        //std::cout << old_sum;
        
        S(num_i, num_j) = S(num_i, num_j)*(-1);
        int new_sum =  S(num_i, num_j)*(S((num_i-1+L)%L, num_j) + S(num_i, (num_j-1+L)%L) + S((num_i+1)%L, num_j) + S(num_i, (num_j+1)%L));
        // std::cout << "\n\n sum_neighbours:";
        // std::cout << new_sum;
        
        S(num_i, num_j) = S(num_i, num_j)*(-1);

        double delta_E = old_sum - new_sum;
       //std::cout << "\n\n delta_E:";
        //std::cout << delta_E;


        double r = disDOUBLE(generator);

        // Check if it accept the change, and flips the spinn
        if(delta_E < 0 || r <= boltzmann_factor[delta_E]){
            S(num_i, num_j) = S(num_i, num_j)*(-1);
            E_tot += delta_E;
        }

    }

    void Ising_model::update(){
            //E = 0;
            //exp_epsilon = 0;
            //E = 0;
            //M = 0;
            E = tot_energy(S);
            M = tot_magnetization(S);
            //exp_m = 0;
            for(int j = 0; j < N; j++){
                MCMC();
                //E_col += E_tot;

                E += tot_energy(S);
                M += tot_magnetization(S);
                }  

            //E = tot_energy(S);
            epsilon = tot_energy(S)/N;
            m = M/(N);
            //E = E/N;
            exp_epsilon = E/(N*N);
            exp_m = M/(N*N);
            // Update expectation values
            //E = E_col;
            //e = E_col/(N);
            //m = abs(B_col)/(N);
            spes_heat = 1./N*1./pow(T,2)*(((E*E)/N)-(pow((E/N),2)));
            suscept = 1./N*1./T*((abs(M*M)/(N*N))-(pow(abs(M)/N, 2)));
    }
    






    



    


