//
// Created by Pushkar Kumar Jain on 3/21/17.
//
#include "Eigen/Dense"
#include "main.h"

Eigen::MatrixXd compute_dj(const Eigen::MatrixXd& xt, double x_g, double y_g, int nnode, int ndim)
{
    Eigen::MatrixXd dfxi(nnode,ndim);
    Eigen::MatrixXd dj(ndim,ndim);
    dfxi = compute_dfxi(x_g,y_g, nnode, ndim);
    dj = xt * dfxi;
    return dj;
}