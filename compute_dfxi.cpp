//
// Created by Pushkar Kumar Jain on 3/21/17.
//
#include "Eigen/Dense"

Eigen::MatrixXd compute_dfxi(double zi, double eta, int nnode, int ndim) {
    Eigen::MatrixXd dfxi(nnode, ndim);
    dfxi(0, 0) = -0.25 * (1 - eta);
    dfxi(1, 0) = 0.25 * (1 - eta);
    dfxi(2, 0) = 0.25 * (1 + eta);
    dfxi(3, 0) = -0.25 * (1 + eta);
    dfxi(0, 1) = -0.25 * (1 - zi);
    dfxi(1, 1) = -0.25 * (1 + zi);
    dfxi(2, 1) = 0.25 * (1 + zi);
    dfxi(3, 1) = 0.25 * (1 - zi);

    return dfxi;
}