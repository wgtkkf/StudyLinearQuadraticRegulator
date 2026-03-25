#include <iostream>    // need this
#include "include/function.h"
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