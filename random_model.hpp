
#include <iostream>
#include <armadillo>


std::mt19937 generator2 (123);
std::uniform_int_distribution<int> disINT2(0.0, 10.);
std::uniform_real_distribution<double> disDOUBLE2(0.0, 10.);

//Declaring function
arma::mat random_model(arma::mat matrix, double L);


arma::mat random_model(arma::mat matrix, double L){

    for(double i = 0; i < L; i++){
        for(double j = 0;j < L; j++){

            double rand_num = disDOUBLE2(generator2);
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
