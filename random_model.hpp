
#include <iostream>
#include <armadillo>


std::mt19937 generator2 (123);
std::uniform_real_distribution<double> disDOUBLE2(0.0, 10.);

//Declaring function
arma::mat random_model(arma::mat matrix, double L);

// Creating a matrix with random spin states
arma::mat random_model(arma::mat matrix, double L){
    // Go trough every position in the matrix
    for(double i = 0; i < L; i++){
        for(double j = 0;j < L; j++){
            
            // Generate random number between 0 and 10
            double rand_num = disDOUBLE2(generator2);

            if(rand_num < 5){
                // Set spin to -1
                matrix(i,j) = -1;
            }
            else{
                // Set spin to 1
                matrix(i,j) = 1;
            }
        }
    }
    return matrix;
}
