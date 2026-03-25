#ifndef _FUNCTION_H_
#define _FUNCTION_H_
#include <Eigen/Dense> /* this is for 2x2 matrix */

struct MatVec {
    Eigen::VectorXd bmat{};
    Eigen::MatrixXd Amat{};    /* Just declare type, NODE_TOTAL is necessary but is defined in the class below */
    Eigen::VectorXd xvector{}; /* Just declare type, NODE_TOTAL is necessary but is defined in the class below */
    Eigen::VectorXd uvector{};
};

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

#endif // _FUNCTION_H_