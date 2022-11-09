
#include <iostream>
#include <armadillo>

//Declaring function
arma::mat random_model(arma::mat matrix, double L);


arma::mat random_model(arma::mat matrix, double L){

    for(double i = 1; i < L+1; i++){
        for(double j = 1;j < L; j++){

            double rand_num = rand() % 10 + 1;
            //std::cout << "\n\n first";
            //std::cout << rand_num;

            if(rand_num < 5){
                //std::cout << "\n\n second";
                //std::cout << rand_num;
                matrix(i,j) = -1;
            }
            else{

                matrix(i,j) = 1;
            }
        }
    }
    return matrix;
}
