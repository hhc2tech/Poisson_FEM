//
// Created by Pushkar Kumar Jain on 3/21/17.
//
#include "Eigen/Dense"

Eigen::MatrixXd compute_shape(double zi, double eta)
{
    Eigen::MatrixXd N(4,1);
    N(0,0) = 0.25*(1-zi)*(1-eta);
    N(1,0) = 0.25*(1+zi)*(1-eta);
    N(2,0) = 0.25*(1+zi)*(1+eta);
    N(3,0) = 0.25*(1-zi)*(1+eta);
    return N;
}
