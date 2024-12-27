/* This is the *.cpp file that defines all of the functions used for compressible D2Q9 lattice Boltzmann (LB) simulations. This includes ./D2Q9_rect.cpp, which contains typedefs and functions for all
 * simulations (compressible and incompressible alike). */


/* Here, we include all of the packages necessary for the code. We also include ./D2Q9_rect.cpp. */
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
#include "./D2Q9_rect.cpp"
#pragma once


/* Here, we compute the macroscopic quantities (density and the velocity fields) at each time step (in the compressible scheme). Following the suggestion in Krüger (Sec. 3.3.3.3, Pg. 69), we unroll
 * the calculations in momentum space and perform each calculation explicitly. I have a good handle on how the loop in momentum space would be input here; I'm just doing this at Krüger's suggestion.
 * These expressions are taken from the aforementioned section in Krüger (specifically, on Pg. 69, Eq. 3.12). For completeness, we also show what the expression in the k loop would look like. Since
 * we compute these expressions first, we calculate them using the distribution function streamed from the *previous* time step. */
void compute_macroscopic_params_compressible(array3D_double prev_dist_fn, array2D_double &rho_simulated_on_lattice, array2D_double &ux_simulated_on_lattice, array2D_double &uy_simulated_on_lattice) {
  int x_size = prev_dist_fn.size();
  int y_size = prev_dist_fn[0].size();

  for (int i = 0; i < x_size; i++) {
    for (int j = 0; j < y_size; j++) {
      rho_simulated_on_lattice[i][j] = prev_dist_fn[i][j][0] + prev_dist_fn[i][j][1] + prev_dist_fn[i][j][2] + prev_dist_fn[i][j][3] + prev_dist_fn[i][j][4]
                                                             + prev_dist_fn[i][j][5] + prev_dist_fn[i][j][6] + prev_dist_fn[i][j][7] + prev_dist_fn[i][j][8];
      ux_simulated_on_lattice[i][j] = (prev_dist_fn[i][j][1] + prev_dist_fn[i][j][5] + prev_dist_fn[i][j][8] - prev_dist_fn[i][j][3] - prev_dist_fn[i][j][6] - prev_dist_fn[i][j][7])
                                                                                                                                                                       /rho_simulated_on_lattice[i][j];
      uy_simulated_on_lattice[i][j] = (prev_dist_fn[i][j][2] + prev_dist_fn[i][j][5] + prev_dist_fn[i][j][6] - prev_dist_fn[i][j][4] - prev_dist_fn[i][j][7] - prev_dist_fn[i][j][8])
                                                                                                                                                                       /rho_simulated_on_lattice[i][j];

      //Note that to use this loop, we need to add four extra arguments to this function: vect1D_double vel_weight = D2Q9_lattice_weights, vect1D_double &c_x_dir = cx_float,
      //vect1D_double &c_y_dir = cy_float, and array3D_double &u_dot_ci.
      //  q_size = vel_weight.size();
      //  sum_rho = 0.0;
      //  for (int k = 0; k < q_num; k++) {
      //    sum_rho += prev_dist_fn[i][j][k];
      //    rho_simulated_on_lattice[i][j] = sum_rho;
      //    ux_simulated_on_lattice[i][j] = (prev_dist_fn[i][j][k] * cx_float[k])/rho_simulated_on_lattice[i][j];
      //    uy_simulated_on_lattice[i][j] = (prev_dist_fn[i][j][k] * cy_float[k])/rho_simulated_on_lattice[i][j];
      //}
    }
  }
}


/* Here, we calculate the distribution functions (in the compressible scheme). Following the suggestion in Krüger Sec. 3.4.7.5 (Pgs. 92-93), I'm "unrolling" the loop over momentum space, and
 * calculating f_colli][j][k = 1] to f_coll[i][j][k = 9] explicitly. According to Krüger, this makes the calculations faster. These also seem to be more accurate for some reason. These equations are
 * given by expanding Eq. 3.54 (Pg. 82) in Krüger (Sec. 3.4.5). These can be further simplified to the equations given in the aforementioned section in Krüger (specifically, on Pg. 93, Eq. 3.65). For
 * completeness, we also show what the "rolled" expression (using the loop in momentum space) would look like.
 *
 * We note that this is used *both* for the initialisation *and* inside each timestep. We use this up front to initialise the equilibrium distribution function, using the pre-initialised x-velocity,
 * y-velocity, and density arrays. (These were initialised when they were declared already.) The speed of sound is a constant in this setup. */

array3D_double compute_f_eq_compressible(array3D_double &f_eq, array2D_double dens, array2D_double x_vel, array2D_double y_vel,
                                                               array2D_double &vel_sq_simulated_on_lattice, const double sound_speed_sq) {

  int x_num = f_eq.size();
  int y_num = f_eq[0].size();

  for (int i = 0; i < x_num; i++) {
    for (int j = 0; j < y_num; j++) {
      vel_sq_simulated_on_lattice[i][j] = std::pow(x_vel[i][j], 2) + std::pow(y_vel[i][j], 2);

      f_eq[i][j][0] = (4.0/9.0) * dens[i][j] * (1.0 - (vel_sq_simulated_on_lattice[i][j]/(2.0 * sound_speed_sq)));
      f_eq[i][j][1] = (1.0/9.0) * dens[i][j] * 
                                    (1.0 + (x_vel[i][j])/sound_speed_sq + (std::pow(x_vel[i][j], 2))/(2.0 * std::pow(sound_speed_sq, 2)) - (vel_sq_simulated_on_lattice[i][j])/(2.0 * sound_speed_sq));
      f_eq[i][j][2] = (1.0/9.0) * dens[i][j] *
                                    (1.0 + (y_vel[i][j])/sound_speed_sq + (std::pow(y_vel[i][j], 2))/(2.0 * std::pow(sound_speed_sq, 2)) - (vel_sq_simulated_on_lattice[i][j])/(2.0 * sound_speed_sq));
      f_eq[i][j][3] = (1.0/9.0) * dens[i][j] *
                                    (1.0 - (x_vel[i][j])/sound_speed_sq + (std::pow(x_vel[i][j], 2))/(2.0 * std::pow(sound_speed_sq, 2)) - (vel_sq_simulated_on_lattice[i][j])/(2.0 * sound_speed_sq));
      f_eq[i][j][4] = (1.0/9.0) * dens[i][j] *
                                    (1.0 - (y_vel[i][j])/sound_speed_sq + (std::pow(y_vel[i][j], 2))/(2.0 * std::pow(sound_speed_sq, 2)) - (vel_sq_simulated_on_lattice[i][j])/(2.0 * sound_speed_sq));
      f_eq[i][j][5] = (1.0/36.0) * dens[i][j] *
                                   (1.0 + (x_vel[i][j] + y_vel[i][j])/sound_speed_sq + (x_vel[i][j] * y_vel[i][j])/(std::pow(sound_speed_sq, 2)) + (vel_sq_simulated_on_lattice[i][j])/sound_speed_sq);
      f_eq[i][j][6] = (1.0/36.0) * dens[i][j] *
                                   (1.0 - (x_vel[i][j] - y_vel[i][j])/sound_speed_sq - (x_vel[i][j] * y_vel[i][j])/(std::pow(sound_speed_sq, 2)) + (vel_sq_simulated_on_lattice[i][j])/sound_speed_sq);
      f_eq[i][j][7] = (1.0/36.0) * dens[i][j] *
                                   (1.0 - (x_vel[i][j] + y_vel[i][j])/sound_speed_sq + (x_vel[i][j] * y_vel[i][j])/(std::pow(sound_speed_sq, 2)) + (vel_sq_simulated_on_lattice[i][j])/sound_speed_sq);
      f_eq[i][j][8] = (1.0/36.0) * dens[i][j] *
                                   (1.0 + (x_vel[i][j] - y_vel[i][j])/sound_speed_sq - (x_vel[i][j] * y_vel[i][j])/(std::pow(sound_speed_sq, 2)) + (vel_sq_simulated_on_lattice[i][j])/sound_speed_sq);

      //Note that to use this loop, we need to add four extra arguments to this function: vect1D_double weight = D2Q9_lattice_weights, vect1D_double &c_x_dir = cx_float,
      //vect1D_double &c_y_dir = cy_float, and array3D_double &u_dot_ci.
      //  q_len = weight.size();
      //  for (int k = 0; k < q_len; k++) {
      //    u_dot_ci[i][j][k] = (x_vel[i][j] * cx_float[k] + y_vel[i][j] * cy_float[k]);
      //    f_eq[i][j][k] = weight[k] * dens[i][j] *
      //                            (1 + u_dot_ci[i][j][k]/sound_speed_sq + std::pow(u_dot_ci[i][j][k], 2)/(2 * std::pow(sound_speed_sq, 2)) - vel_sq_simulated_on_lattice[i][j]/(2 * sound_speed_sq));
      //  }
    }
  }

  return f_eq;
}