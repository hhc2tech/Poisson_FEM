//
// Created by Pushkar Kumar Jain on 3/21/17.
//
#include "Eigen/Dense"

#ifndef UNTITLED_MAIN_H
#define UNTITLED_MAIN_H

#endif //UNTITLED_MAIN_H

Eigen::MatrixXd compute_shape(double zi, double eta);

Eigen::MatrixXd compute_dfxi(double zi, double eta, int nnode, int ndim);

Eigen::MatrixXd compute_dj(const Eigen::MatrixXd& xt, double x_g, double y_g, int nnode, int ndim);

Eigen::MatrixXd compute_dji(const Eigen::MatrixXd&  xt, double x_g, double y_g, int nnode, int ndim);

Eigen::MatrixXd compute_dfx(const Eigen::MatrixXd& xt, double x_g, double y_g, int nnode, int ndim);

struct coordinate
{
    double x;
    double y;
};

struct connectivity
{
    int node[4];
};
