//
// Created by Pushkar Kumar Jain on 3/21/17.
//
#include "Eigen/Dense"

struct coordinate
{
    double x;
    double y;
};

struct connectivity
{
    int node[4];
};

std::vector<coordinate> set_mesh(double xstart, double ystart, int num_nodes_x, int num_nodes_y, int nj, double h);

std::vector<int> set_bc(int num_nodes_x, int num_nodes_y);

std::vector<connectivity> set_connectivity(int nelem, int num_nodes_x);

Eigen::MatrixXd compute_shape(double zi, double eta);

Eigen::MatrixXd compute_dfxi(double zi, double eta, int nnode, int ndim);

Eigen::MatrixXd compute_dj(const Eigen::MatrixXd& xt, double x_g, double y_g, int nnode, int ndim);

Eigen::MatrixXd compute_dji(const Eigen::MatrixXd&  xt, double x_g, double y_g, int nnode, int ndim);

Eigen::MatrixXd compute_dfx(const Eigen::MatrixXd& xt, double x_g, double y_g, int nnode, int ndim);


