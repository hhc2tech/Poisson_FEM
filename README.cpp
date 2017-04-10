/*! \mainpage Mainpage
 *
 * \section intro_sec Introduction
 *
 * This code uses the Finite element method to solve the Poisson heat equation with a source term
 * on a rectangular domain.
 *
 * The PDE is
 * \f[
 * \Delta u = \sin(\pi x)\ on\ \Omega,
 * u = 0\ on\ \partial\Omega
 * \f]
 *
 * where \f$ \Omega = {(x,y):0<x<2,0<y<1}\f$ and \f$\partial\Omega\f$ is the boundary of \f$\Omega\f$.
 *
 * Bilinear rectangular element is used.
 *
 * \section install_sec Weak form
 *
 * Find \f$u\in H_{0}^{1}\f$ such that
 * \f[
 * \int_{\Omega}\nabla v\cdot\nabla u\ dx=\int_{\Omega}v\sin\left(\pi x\right)dx\ \forall v\in H_{0}^{1}
 * \f]
 *
 * \section steps Steps
 *
 *   -# Set the computational domain - Mesh, connectivity and boundary conditions
 *   -# Set quadrautre integration method
 *   -# Loop over every element to calculate local force and local stiffness
 *   -# Assemble local stiffness and local force to global stiffness \f$(K)\f$ and global force \f$(F)\f$ respectively
 *   -# Apply boundary conditions
 *   -# Solve \f$Kx=F\f$
 */