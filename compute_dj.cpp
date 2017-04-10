/**
 *  @file    compute_dj.cpp
 *  @author  Pushkar Kumar Jain (pushkarjain1991@utexas.edu)
 *  @date    4/9/2017
 *  @version 1.0
 *
 *  @brief Computes Jacobian
 *
 *  @section DESCRIPTION
 *
 *  This is a little program that reads a list of names from
 *  a specified file or from standard input and then sorts
 *  the names in ascending order and prints them to standard
 *  output.
 *
 *  Command line arguments are used to specify where the
 *  list of names should be read from.  If the program
 *  doesn't receive any command line arguments then it reads
 *  the names from standard input. If the program receives
 *  a single command line argument then it reads the names
 *  from the corresponding file.  If more than one command
 *  line argument is specified the program prints a usage
 *  message and terminates.
 *
 */
#include "Eigen/Dense"
#include "main.h"

/**
*   @brief Computes the Jacobian
*
*   @param  xt is an initialized integer variable
*   @param  x_g is an initialized integer variable
*   @param  y_g is an initialized integer variable
*   @param  nnode is an initialized integer variable
*   @param  ndim is an initialized integer variable
*   @return Eigen::MatrixXd
*/
Eigen::MatrixXd compute_dj(const Eigen::MatrixXd &xt, double x_g, double y_g, int nnode, int ndim) {

    Eigen::MatrixXd dfxi(nnode, ndim);
    Eigen::MatrixXd dj(ndim, ndim);
    dfxi = compute_dfxi(x_g, y_g, nnode, ndim);
    dj = xt * dfxi;
    return dj;
}