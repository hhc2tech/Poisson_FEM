

/**
 *  @file    main.cpp
 *  @author  Pushkar Kumar Jain (pushkarjain1991@utexas.edu)
 *  @date    4/9/2017
 *  @version 1.0
 *
 *  @brief The main program to feed inputs and run
 *
 */

#include <iostream>
#include <math.h>
#include <vector>
#include "Eigen/Dense"
#include <unsupported/Eigen/MatrixFunctions>
#include "main.h"


int main() {
    double zi_1[2], eta_1[2];
    double J_dx_dy;
    int IG, JG;
    int iter1, iter2;

    // Input conditions for rectangular domain
    const double x_start = 0.0; // Initial point in domain in x direction
    const double x_end = 2.0; // Final point in domain in x direction
    const double y_start = 0.0; // Initial point in domain in y direction
    const double y_end = 1.0; // Final point in domain in y direction
    const double h = 0.25; // Domain spacing

    // FEM input conditions
    // Choosing rectangular elements parameters
    const int nnode = 4; // Rectangular elements have 4 nodes
    const int ndim = 2; // 2 dimensions x and y

    //Integration parameters
    zi_1[0] = -1.0 / sqrt(3.0);
    zi_1[1] = 1.0 / sqrt(3.0);
    eta_1[0] = -1.0 / sqrt(3.0);
    eta_1[1] = 1.0 / sqrt(3.0);
    const int num_gauss_quad_points = 4; // Number of Gauss points for integration;

    std::vector<coordinate> gauss_points(num_gauss_quad_points);
    int counter = 0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            gauss_points[counter].x = zi_1[i];
            gauss_points[counter].y = eta_1[j];
            counter++;
        }
    }

    const int num_nodes_x = static_cast<int> (floor((x_end - x_start) / h) + 1); // Number of nodes in x
    const int num_nodes_y = static_cast<int> (floor((y_end - y_start) / h) + 1); // Number of nodes in y
    const int nj = num_nodes_x * num_nodes_y; // Total number of nodes
    const int nelem = (num_nodes_x - 1) * (num_nodes_y - 1); // Total number of elements

    // Setting up the mesh
    std::vector<coordinate> mesh(nj);
    mesh = set_mesh(x_start, y_start, num_nodes_x, num_nodes_y, nj, h);

    // Setting boundary conditions
    std::vector<int> id(nj);
    id = set_bc(num_nodes_x, num_nodes_y);

    // Setting connectivity
    std::vector<connectivity> inode(nelem);
    inode = set_connectivity(nelem, num_nodes_x);

    Eigen::MatrixXd global_k(nj, nj); // Initializing global stiffness
    global_k.setZero();

    Eigen::VectorXd global_f(nj); // Initializing global force
    global_f.setZero();

    for (int element_num = 0; element_num < nelem; ++element_num) {
        Eigen::MatrixXd local_k(nnode, nnode); // Initializing local force
        local_k.setZero();

        Eigen::VectorXd local_f(nnode); // Initializing local force
        local_f.setZero();

        Eigen::MatrixXd xt(ndim, nnode); // xt is ndim X nnode matrix with x,y coordinates of an element
        for (int i = 0; i < nnode; ++i) {
            xt(0, i) = mesh[inode[element_num].node[i]].x;
            xt(1, i) = mesh[inode[element_num].node[i]].y;
        }

        // Integrating with gaussian quadrature
        for (int i = 0; i < num_gauss_quad_points; ++i) {
            Eigen::MatrixXd grad_N(nnode, ndim);
            Eigen::MatrixXd N(nnode, 1);
            Eigen::MatrixXd J(ndim, ndim);

            grad_N = compute_dfx(xt, gauss_points[i].x, gauss_points[i].y, nnode, ndim);
            N = compute_shape(gauss_points[i].x, gauss_points[i].y);
            J = compute_dj(xt, gauss_points[i].x, gauss_points[i].y, nnode, ndim);
            J.transposeInPlace();
            J_dx_dy = J.determinant();

            for (int j = 0; j < ndim; ++j) {
                local_k = local_k + J_dx_dy * grad_N.col(j) * grad_N.col(j).transpose();
            }

            local_f = local_f + sin(xt.row(0) * N.col(0)) * J_dx_dy * N;
        }

        // Assembly to global K and F
        for (int i = 0; i < nnode; ++i) {
            IG = inode[element_num].node[i];
            for (int j = 0; j < nnode; ++j) {
                JG = inode[element_num].node[j];
                global_k(IG, JG) = global_k(IG, JG) + local_k(i, j);
            }
            global_f(IG) = global_f(IG) + local_f(i);
        }
    }

    // Apply BC
    for (int eq_num = 0; eq_num < nj; ++eq_num) {
        if (id[eq_num] != 0) {
            global_k.row(eq_num).setZero();
            global_k.col(eq_num).setZero();
            global_k(eq_num, eq_num) = 1.0;
            global_f(eq_num) = 0.0;
        }
    }

    // Solving system of equations
    Eigen::FullPivLU<Eigen::MatrixXd> LU(global_k);
    Eigen::VectorXd solution = LU.solve(global_f);

    std::cout << "Solution" << std::endl;
    std::cout << solution;

    return 0;
}
