#include <iostream>    // need this
#include "include/function.h"
#include "include/datatypes.h"
#include <cmath>       // for mathemetical operator such as power function
#include <format>      // for printing out
#include <iomanip>     // For formatting
#include <Eigen/Dense> // For partialPivLu


/* Message methods */
void Message::start()
{
    std::cout << "### Calculation started.   ###\n";
}

void Message::end()
{
    std::cout << "### Calculation completed. ###\n";
}

void Display::show(Eigen::MatrixXd* Amat, Eigen::VectorXd* bmat)
{
    std::cout << "--- Matrix Amat ---" << std::endl;
    std::cout << *Amat << std::endl;
    std::cout << "-------------------" << std::endl;

    std::cout << "--- Matrix Bmat ---" << std::endl;
    std::cout << *bmat << std::endl;
    std::cout << "-------------------" << std::endl;
}

/* This is a constructor */
MatVecGeneration::MatVecGeneration(){
    inputs.Amat.setZero();
    inputs.Bmat.setZero();
    inputs.Qmat.setZero();
    inputs.Rmat.setZero();
    inputs.Pmat.setZero();
}

double MatVecGeneration::DefineInputs(){    
    inputs.Amat(0, 1) = 1.0;    

    return 0;
}