#include <iostream>
#include <math.h>
#include <vector>
#include "Eigen/Dense"
#include <unsupported/Eigen/MatrixFunctions>
#include "main.h"

int main() {
    //const double pi = std::acos(-1.0);
    double x_start, x_end, y_start, y_end, h;
    int nnode, ndim;
    int num_gauss_quad_points;
    double zi_1[2], eta_1[2];
    int num_nodes_x, num_nodes_y, nj, nelem;
    double x_coordinate, y_coordinate;
    double J_dx_dy;
    int IG, JG;
    int iter1, iter2;

    // Input conditions for rectangular domain
    x_start = 0.0;
    x_end = 2.0;
    y_start = 0.0;
    y_end = 1.0;
    h = 0.25;

    // FEM input conditions
    // Choosing rectangular elements
    nnode = 4;
    ndim = 2;

    //Integration input conditions
    zi_1[0] = -1.0 / sqrt(3.0);
    zi_1[1] = 1.0 / sqrt(3.0);
    eta_1[0] = -1.0 / sqrt(3.0);
    eta_1[1] = 1.0 / sqrt(3.0);
    num_gauss_quad_points = 4;


    std::vector<coordinate> gauss_points(num_gauss_quad_points);
    int counter = 0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            gauss_points[counter].x = zi_1[i];
            gauss_points[counter].y = eta_1[j];
            counter++;
        }
    }

    num_nodes_x = static_cast<int> (floor((x_end - x_start) / h) + 1);
    num_nodes_y = static_cast<int> (floor((y_end - y_start) / h) + 1);
    nj = num_nodes_x * num_nodes_y;
    nelem = (num_nodes_x - 1) * (num_nodes_y - 1);


    std::vector<coordinate> mesh(nj);
    y_coordinate = y_start;
    counter = 0;
    for (int i = 0; i < num_nodes_y; ++i) {
        x_coordinate = x_start;

        for (int j = 0; j < num_nodes_x; ++j) {
            mesh[counter].y = y_coordinate;
            mesh[counter].x = x_coordinate;
            x_coordinate += h;
            counter++;
        }
        y_coordinate += h;
    }

    // Assign BC
    std::vector<int> id(nj);
    counter = 0;
    for (int i = 0; i < num_nodes_y; ++i) {
        for (int j = 0; j < num_nodes_x; ++j) {
            (i == 0 || i == num_nodes_y - 1 || j == 0 || j == num_nodes_x - 1) ? id[counter] = 1 : id[counter] = 0;
            counter++;
        }
    }

    // Setting the connectivity
    std::vector<connectivity> inode(nelem);
    iter1 = 0;
    iter2 = 0;
    for (int element_num = 0; element_num < nelem; ++element_num) {
        inode[element_num].node[0] = iter1 + iter2 * (num_nodes_x);
        inode[element_num].node[1] = inode[element_num].node[0] + 1;
        inode[element_num].node[2] = inode[element_num].node[1] + num_nodes_x;
        inode[element_num].node[3] = inode[element_num].node[2] - 1;
        iter1++;
        if ((element_num + 1) % (num_nodes_x - 1) == 0) {
            iter1 = 0;
            iter2++;
        }
    }

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
