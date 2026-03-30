#ifndef _FUNCTION_H_
#define _FUNCTION_H_
#include <Eigen/Dense> /* this is for 2x2 matrix */
#include "datatypes.h" /* the same directory, include struct type */

class Message
{
    public:
        void start();
        void end();
};

class Display /* Finite Element class */
{
    public:                        
        void show(Eigen::MatrixXd* Amat, Eigen::VectorXd* bmat);
};

class MatVecGeneration{
    public:
        MatVec inputs;
        MatVecGeneration();     /* This is a constructor */
        double DefineInputs();

};

#endif // _FUNCTION_H_