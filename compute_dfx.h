//
// Created by Pushkar Kumar Jain on 3/21/17.
//
#include "Eigen/Dense"
#include "main.h"

#ifndef UNTITLED_COMPUTE_DFX_H
#define UNTITLED_COMPUTE_DFX_H

#endif //UNTITLED_COMPUTE_DFX_H

Eigen::MatrixXd compute_dfxi(double zi, double eta, int nnode, int ndim);
Eigen::MatrixXd compute_dji(const Eigen::MatrixXd&  xt, double x_g, double y_g, int nnode, int ndim);