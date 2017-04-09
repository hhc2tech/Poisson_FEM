//
// Created by Pushkar Kumar Jain on 3/21/17.
//
#include "Eigen/Dense"
#include "main.h"


Eigen::MatrixXd compute_dji(const Eigen::MatrixXd&  xt, double x_g, double y_g, int nnode, int ndim)
{
    double denom;
    Eigen::MatrixXd dji(ndim,ndim);
    Eigen::MatrixXd temp_mat(ndim,ndim);
    temp_mat = compute_dj(xt, x_g, y_g, nnode, ndim);
    dji(0,0) = temp_mat(1,1);
    dji(1,1) = temp_mat(0,0);
    dji(1,0) = -temp_mat(1,0);
    dji(0,1) = -temp_mat(0,1);
    denom = 1.0/temp_mat.determinant();
    dji = denom*dji;
    return dji;
}