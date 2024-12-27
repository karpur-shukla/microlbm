/* This is the *.cpp file that defines all of the functions used for incompressible D2Q9 lattice Boltzmann (LB) simulations. This includes ./D2Q9_rect.cpp, which contains typedefs and functions for
 * all simulations (compressible and incompressible alike). */


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


/* Here, we compute the macroscopic quantities (density and the velocity fields) at each time step (in the incompressible scheme). Following the suggestion in Krüger (Sec. 3.3.3.3, Pg. 69), we unroll
 * the calculations in momentum space and perform each calculation explicitly. I have a good handle on how the loop in momentum space would be input here; I'm just doing this at Krüger's suggestion.
 * These expressions are taken from the aforementioned section in Krüger (specifically, on Pg. 69, Eq. 3.12). For completeness, we also show what the expression in the k loop would look like. Since
 * we compute these expressions first, we calculate them using the distribution function streamed from the *previous* time step. */
void compute_macroscopic_params_incompressible(array3D_double prev_dist_fn, array2D_double &rho, array2D_double &ux, array2D_double &uy) {
  int x_size = prev_dist_fn.size();
  int y_size = prev_dist_fn[0].size();

  for (int i = 0; i < x_size; i++) {
    for (int j = 0; j < y_size; j++) {
      rho[i][j] = prev_dist_fn[i][j][0] + prev_dist_fn[i][j][1] + prev_dist_fn[i][j][2] + prev_dist_fn[i][j][3] + prev_dist_fn[i][j][4]
                                                                + prev_dist_fn[i][j][5] + prev_dist_fn[i][j][6] + prev_dist_fn[i][j][7] + prev_dist_fn[i][j][8];
      ux[i][j] = prev_dist_fn[i][j][1] + prev_dist_fn[i][j][5] + prev_dist_fn[i][j][8] - prev_dist_fn[i][j][3] - prev_dist_fn[i][j][6] - prev_dist_fn[i][j][7];
      uy[i][j] = prev_dist_fn[i][j][2] + prev_dist_fn[i][j][5] + prev_dist_fn[i][j][6] - prev_dist_fn[i][j][4] - prev_dist_fn[i][j][7] - prev_dist_fn[i][j][8];

      //Note that to use this loop, we need to add four extra arguments to this function: vect1D_double vel_weight = D2Q9_lattice_weights, vect1D_double &c_x_dir = cx_float,
      //vect1D_double &c_y_dir = cy_float, and array3D_double &u_dot_ci.
      //  q_size = vel_weight.size();
      //  sum_rho = 0.0;
      //  for (int k = 0; k < q_num; k++) {
      //    sum_rho += prev_dist_fn[i][j][k];
      //    rho_simulated_on_lattice[i][j] = sum_rho;
      //    ux_simulated_on_lattice[i][j] = prev_dist_fn[i][j][k] * cx_float[k];
      //    uy_simulated_on_lattice[i][j] = prev_dist_fn[i][j][k] * cy_float[k];
      //}
    }
  }
}


/* Here, we calculate the distribution functions (in the incompressible scheme). Unlike everything else, I haven't yet unrolled these expressions, due to a lack of time. This uses Eq. 5.35 on Pg. 192
 * of Krüger.
 *
 * We note that this is used *both* for the initialisation *and* inside each timestep. We use this up front to initialise the equilibrium distribution function, using the pre-initialised x-velocity,
 * y-velocity, and density arrays. (These were initialised when they were declared already.) The speed of sound is a constant in this setup. */
void compute_f_eq_incompressible(array2D_double dens, array2D_double x_vel, array2D_double y_vel, vect1D_double weight, vect1D_double c_in_x_dir, vect1D_double c_in_y_dir,
                                                      array3D_double &f_eq, array3D_double &u_dot_ci, const double sound_speed_sq) {
  int x_num = f_eq.size();
  int y_num = f_eq[0].size();
  int q_len = weight.size();
  for (int i = 0; i < x_num; i++) {
    for (int j = 0; j < y_num; j++) {
      for (int k = 0; k < q_len; k++) {
        u_dot_ci[i][j][k] = (x_vel[i][j] * c_in_x_dir[k] + y_vel[i][j] * c_in_y_dir[k]);
        f_eq[i][j][k] = weight[k] * (dens[i][j] + u_dot_ci[i][j][k]/ sound_speed_sq);
      }
    }
  }
}