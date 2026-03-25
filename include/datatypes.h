#pragma once
#include <Eigen/Dense> /* this is for 2x2 matrix */

struct MatVec {
    static constexpr int size4 = 4;
    static constexpr int size1 = 1;

    Eigen::MatrixXd Amat{size4, size4};
    Eigen::MatrixXd Bmat{size4, size1};
    Eigen::MatrixXd Qmat{size4, size4};
    Eigen::MatrixXd Rmat{size1, size1};
    Eigen::MatrixXd Pmat{size4, size4};
};
