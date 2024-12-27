/* This is the *.cpp file that defines all of the functions used in every square-coordinate D2Q9 lattice Boltzmann (LB) simulation, except for the time steppers. At present, this is primarily used to
 * simulate fully-developed 2D Couette flow, fully-developed 2D Poiseuille flow, and the 2D lid-driven cavity problem. All of these simulations are isothermal. 
 * 
 * Here, we primarily reference Krüger et al., The Lattice Boltzmann Method, Springer Nature, Cham, Switzerland (2017) and D. Wilde et al., Int. J. Numer. Meth. Fluids 90, 156 (2019),
 * DOI:10.1002/fld.4716. Other references will be indicated when applicable. 
 * 
 * At present, we use the following setup:
 *   Collision operator:    BGK approximation
 *   Lattice structure:     D2Q9
 *   1st order BCs:         Periodic in x; wet-node bounce-back (c.f. Krüger Pg. 177) in y at walls
 *   2nd order (noneq) BCs: Zou-He (nonequilibrium bounce back) (c.f. Krüger Pgs. 196-199). The equations for the nonequilibrium extrapolation method (c.f. Krüger Pgs. 194-195) are also included for
 *                          reference.
 *   Relaxation scheme:     Single relaxation time
 *   Added forces:          None
 *   Number of phases:      1
 * 
 * Planned extensions are:
 * - Extension of the collision operator to RK4 and TVD-RK3 (as in A. Shirasat, S. Nayak, and D. Patil, Phys. Rev. E 106, 025314 (2022), DOI:10.1103/PhysRevE.106.025314).
 * - Extension of relaxation time structure to Shakhov and/or multiple relaxation times.
 * - Multiphase (via Shan-Chen as a first step).
 * - Use of immersed-boundary methods to simulate obstacles.
 * - Extension of the spatial grid discretisation to finite-volume problems. 
 * - Reading inputs from a file.
 * - Creating a MAKEFILE. */


/* Here, we include all of the packages necessary for the code. */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <numeric>
#include "D2Q9_rect.h"
#pragma once


/* Here, we define the 2D array class and 3D array class as specific types. This is to simplify the expressions when we initialise the 1-particle probability distribution functions, which are 3D
 * arrays. With this, we can avoid writing std::vectors inside other std::vectors, which quickly becomes impossible to read or deal with. For reference, before I did the typedef, one of these looked
 * like: std::vector<std::vector<std::vector<double>>> f(x_len, std::vector<std::vector<double>>(y_len, std::vector<double>(q_num))). */
typedef std::vector<double> vect1D_double;
typedef std::vector<std::vector<double>> array2D_double;
typedef std::vector<std::vector<std::vector<double>>> array3D_double;


/* These are the simulation parameters, including some "extra" simulation parameters. Krüger Sec. 3.3.3.4 (Pgs. 69-70) suggests to also define cs⁴ (speed of sound to the fourth) up front. Here, we
 * don't use that, so I'm just defining cs². */
const int x_len = 2;                           // number of x grid points. System is periodic in x.
const int y_len = 21;                           // number of y grid points
const int t_steps = 10001;                   // number of time steps. Since I update the densities and velocities at the beginning of the simulation, I'm doing one extra timestep.
const double BGK_relax_time_current_FECN = 1.0; // BGK relaxation time for forward Euler / Crank-Nicolson
const double BGK_relax_time_current_RK4 = 1.0;  // BGK relaxation time for 4th-order Runge-Kutta
const double BGK_relax_time_current_AM3 = 1.0;  // BGK relaxation time for 3rd-order Adams-Moulton
const double BGK_relax_time_current_BDF2 = 1.0; // BGK relaxation time for 2nd-order backwards differentiation
const double rho_global_start = 1.0;            // initial density
const double cs_sq = 1.0/3.0;                   // speed of sound squared in lattice units
const double ux_init_global = 0.0;              // global initial velocity in the x-direction
const double uy_init_global = 0.0;              // global initial velocity in the y-direction


/* Here, we define quantities that will be used for Couette and Poiseuille flow.  */
double u_max_centre = 0.0;                  // centre-line velocity (used for pressure-driven Poiseuille flow)
double u_wall_top = 0.0;                    // top wall velocity (used for Couette flow)
double u_wall_bottom = 0.0;                 // bottom wall velocity (used for Couette flow)
const double u_wall_top_reference = 1.0;    // reference top wall velocity (used for Couette flow over an order of magnitude)
const double u_wall_bottom_reference = 0.0; // reference top wall velocity (used for Couette flow over an order of magnitude)
const double u_max_centre_reference = 1.0;  // reference centre-line velocity (used for Poiseuille flow)


/* These are "derived" parameters. For the initialisation, I'm defininig vel_init_sq (the square of the *total* initial velocity; i.e. (ux_init_global)² + (uy_init_global)²). Defining vel_init_sq
 * just lets the initialisation be more comprehensibleThe expression for the viscosity used here is given by kinematic_visc = rho * cs_sq (BGK_relax_time_current - dt/2), where dt is the timestep.
 * This is Eq. 4.17 in Sec. 4.1.4 of Krüger (Pg. 112); this is a generalisation of Eq. 3.5 in Sec. 3.2.1 of Krüger (Pg. 65), which doesn't include rho. However, we note that an alternate expression
 * for the (nondimensionalised) viscosity is given by Eq. 7.14 in Sec. 7.2.1.1 of Krüger (Pg. 273) as kinematic_visc = cs_sq * (BGK_relaxation_time - 1/2) * (dx)²/dt. 
 * 
 * The expressions for the BGK relaxation times for the "explicitised" versions of the 3rd-order accurate Adams-Moulton method and the 2nd-order accurate backwards differentiation formula method are
 * given in D. Wilde et al., Pg. 159 (the 4th page of the article). */
const double vel_init_sq = std::pow(ux_init_global, 2) + std::pow(uy_init_global, 2);                  // square of the magnitude of the initial velocity
const double BGK_relax_time_prev_AM3 = 13.0 * BGK_relax_time_current_AM3;                              // BGK relaxation time for previous timestep for 3rd-order Adams-Moulton
const double BGK_relax_time_prev_BDF2 = 4.0 * BGK_relax_time_current_BDF2;                             // BGK relaxation time for previous timestep for 2nd-order backwards differentiation

const double inverse_relaxation_time_FECN = 1.0/BGK_relax_time_current_FECN;                           // BGK inverse relaxation time for forward Euler / Crank-Nicolson (lowercase omega, NOT capital omega)
const double inverse_relaxation_time_AM3 = 1.0/BGK_relax_time_current_AM3;                             // BGK inverse relaxation time for 3rd-order Adams-Moulton (lowercase omega, NOT capital omega)
const double inverse_relaxation_time_BDF2 = 1.0/BGK_relax_time_current_BDF2;                           // BGK inverse relaxation time for 2nd-order backwards differentiation (lowercase omega, NOT capital omega)
const double inverse_relaxation_time_RK4 = 1.0/BGK_relax_time_current_RK4;                             // BGK inverse relaxation time for 4th-order Runge-Kutta (lowercase omega, NOT capital omega)

const double kinemat_visc_FECN = cs_sq * rho_global_start * (BGK_relax_time_current_FECN - 1.0 / 2.0); // kinematic viscosity for forward Euler / Crank-Nicolson
const double kinemat_visc_AM3 = cs_sq * rho_global_start * (BGK_relax_time_current_AM3 - 1.0/2.0);     // kinematic viscosity for 3rd-order Adams-Moulton
const double kinemat_visc_BDF2 = cs_sq * rho_global_start * (BGK_relax_time_current_BDF2 - 1.0/2.0);   // kinematic viscosity for 2nd-order backwards differentiation
const double kinemat_visc_RK4 = cs_sq * rho_global_start * (BGK_relax_time_current_RK4 - 1.0/2.0);     // kinematic viscosity for 2nd-order backwards differentiation


/* These are the lattice parameters. We define the x-direction and y-direction lattice speeds as both a float vector and an int vector. The float vector goes into equations in which the lattice
 * projections (the c_i's, i.e., the x-direction and y-direction lattice speeds) are used. (Usually, this involves expressions involving the dot product of u and c.) */
const int q_num = 9; // number of velocity directions (Q in DnQm). Here, Q = 9, with q = 0 as the self-velocity.

const std::vector<double> D2Q9_lattice_weights = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0}; // lattice weights

const std::vector<int> cx_int = {0, 1, 0, -1,  0, 1, -1, -1,  1};                        // x-direction lattice speeds
const std::vector<int> cy_int = {0, 0, 1,  0, -1, 1,  1, -1, -1};                        // y-direction lattice speeds
const std::vector<double> cx_float = {0.0, 1.0, 0.0, -1.0,  0.0, 1.0, -1.0, -1.0,  1.0}; // x-direction lattice speeds
const std::vector<double> cy_float = {0.0, 0.0, 1.0,  0.0, -1.0, 1.0,  1.0, -1.0, -1.0}; // y-direction lattice speeds

/* The lattice velocity directions above follow Krüger Pg. 86. These are defined by:
 *
 * index:   0  1  2  3  4  5  6  7  8
 * ----------------------------------
 * x:       0 +1  0 -1  0 +1 -1 -1 +1
 * y:       0  0 +1  0 -1 +1 +1 -1 -1
 *
 * 6 2 5
 *  \|/    ↑ y
 * 3-0-1   |
 *  /|\    --→ x
 * 7 4 8 */


 /* This defines the distribution function arrays. Since the distribution functions live in phase space, we have one f for each value of x, y, AND q; i.e., at each node, we have q_num different
  * distribution functions. (Here, q_num = 9.) Thus, if we want to consider each f, f_eq, and f_prop (the distribution function, equilibrium distribution function, and propagated distribution
  * function respectively), we need to make three rank-3 arrays. */
array3D_double f_coll(x_len, array2D_double(y_len, vect1D_double(q_num, 0.0)));      // this is f_i in Krüger Ch. 4; f^[+1] in Wilde et al.
array3D_double f_eq(x_len, array2D_double(y_len, vect1D_double(q_num, 0.0)));        // this is f_i^eq in Krüger Ch. 4; f^[0]_eq in Wilde et al.
array3D_double f_prop(x_len, array2D_double(y_len, vect1D_double(q_num, 0.0)));      // this is f_i^* in Krüger Ch. 4; f^[0] in Wilde et al.
array3D_double f_prop_prev(x_len, array2D_double(y_len, vect1D_double(q_num, 0.0))); // this is f_i^* in Krüger Ch. 4; f^[-1] in Wilde et al.
array3D_double f_eq_prev(x_len, array2D_double(y_len, vect1D_double(q_num, 0.0)));   // this is f_i^eq in Krüger Ch. 4; f^[-1]_eq in Wilde et al.


/* Here, we define the macroscopic quantities (for this problem, the density and velocity fields). In their declaration, we initialise them up-front with their initial values (rho_global_start for
 * the initial density, and ux_init_global and uy_init_global respectively for the initial x- and y-velocities). This will allow us to immediately initialise the equilibrium distribution function
 * f_eq using the exact same compute_f_eq function that we use to compute f_eq at each timestep. */
array2D_double rho_simulated_on_lattice(x_len, vect1D_double(y_len, rho_global_start)); // density array
array2D_double ux_simulated_on_lattice(x_len, vect1D_double(y_len, ux_init_global));    // x-direction velocity array
array2D_double uy_simulated_on_lattice(x_len, vect1D_double(y_len, uy_init_global));    // y-direction velocity array
array2D_double vel_sq_simulated_on_lattice(x_len, vect1D_double(y_len, vel_init_sq));   // array to store the magnitude of the macroscopic velocity squared
array3D_double u_dot_ci(x_len, array2D_double(y_len, vect1D_double(q_num, 0.0)));       // array to store the inner product of u (the macroscopic velocity) and c (the lattice velocity projections)


/* Here, we perform the streaming step, without applying the boundary condition at the walls, but with applying periodicity. The boundary condition at the wall shifts everything over, and this will
 * be applied in a different function. We note that this is fine, since we can impose the boundary conditions at a different step. A discussion of applying periodic boundary conditions everywhere is
 * given at the Palabos forum post at https://palabos-forum.unige.ch/t/circshift-in-cavity-m/685 (mirrored at
 * https://web.archive.org/web/20230402015817/https://palabos-forum.unige.ch/t/circshift-in-cavity-m/685). */
array3D_double streaming_periodic_x_and_y(std::vector<int> cx_vec, std::vector<int> cy_vec, array3D_double &f_prop, array3D_double current_coll_dist_fn) {
  int new_x;
  int new_y;

  int x_range = f_prop.size();
  int y_range = f_prop[0].size();
  int q_range = f_prop[0][0].size();

  for (int i = 0; i < x_range; i++) {
    for (int j = 0; j < y_range; j++) {
      for (int k = 0; k < q_range; k++) {
        new_x = (i + cx_vec[k] + x_range) % x_range;
        new_y = (j + cy_vec[k] + y_range) % y_range;
        f_prop[new_x][new_y][k] = current_coll_dist_fn[i][j][k];
      }
    }
  }

  return f_prop;
}


/* Here, we reset the values of the macroscopic rho, x-velocity, and y-velocity, in order to start each simulation from a "clean" set of arrays. */
void reset_macroscopic_parameters(array2D_double density_from_previous_sim, array2D_double ux_from_previous_sim, array2D_double uy_from_previous_sim, array2D_double vel_sq_from_previous_sim,
                                                                            const double initial_density, const double initial_global_x_vel, const double initial_global_y_vel) {
  
  const double initial_global_vel_sq = std::pow(initial_global_x_vel, 2) + std::pow(initial_global_x_vel, 2);

  int sim_size_x_dir = density_from_previous_sim.size();
  int sim_size_y_dir = density_from_previous_sim[0].size();

  for (int i = 0; i < sim_size_x_dir; i++) {
    for (int j = 0; j < sim_size_y_dir; j++) {
      density_from_previous_sim[i][j] = initial_density;
      ux_from_previous_sim[i][j] = initial_global_x_vel;
      uy_from_previous_sim[i][j] = initial_global_y_vel;
      vel_sq_from_previous_sim[i][j] = initial_global_vel_sq;
    }
  }
}


/* Here, we define the function that prints a rank-2 matrix to screen. */
void print_2D_array(array2D_double input_mat) {
  int num_cols = input_mat.size();
  int num_rows = input_mat[0].size();

  for (int i = 0; i < num_cols; i++) {
    std::cout << "[";
    for (int j = 0; j < num_rows - 1; j++) {
      std::cout << input_mat[i][j] << " ";
    }
    std::cout << input_mat[i][num_rows - 1] << "]" << std::endl;
  }
}


/* Here, we define the function that saves a rank-2 matrix to an external *.CSV file, with a given input filename. */
void save_2D_array_to_csv(array2D_double input_matrix, std::string path, std::string filename) {
  int qty_cols = input_matrix.size();
  int qty_rows = input_matrix[0].size();

  std::string full_filename = path + filename + ".csv";

  std::ofstream saved_array(full_filename);

  for (int i = 0; i < qty_cols; i++) {
    for (int j = 0; j < qty_rows - 1; j++) {
      saved_array << input_matrix[i][j] << ",";
    }
    saved_array << input_matrix[i][qty_rows - 1] << std::endl;
  }
}