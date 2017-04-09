//
// Created by Pushkar Kumar Jain on 3/21/17.
//
#include "Eigen/Dense"
#include "main.h"

Eigen::MatrixXd compute_dfx(const Eigen::MatrixXd& xt, double x_g, double y_g, int nnode, int ndim)
{
    Eigen::MatrixXd dfx(nnode,ndim);
    dfx=compute_dfxi(x_g, y_g, nnode, ndim)*compute_dji(xt, x_g, y_g, nnode, ndim);
    return dfx;
}
