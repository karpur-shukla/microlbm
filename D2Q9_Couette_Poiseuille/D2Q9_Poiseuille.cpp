/* This is the *.cpp file that defines all of the functions used to simulate 2D rectangular Poiseuille flow using the D2Q9 lattice Boltzmann (LB) scheme. This includes ./D2Q9_rect.cpp, which contains
 * typedefs and functions for all simulations (compressible and incompressible alike). Here, we note that the Poiseuille flow problem involves a 'halo' region at x = 0 and x = x_len - 1. */


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
#include "./../D2Q9_general/D2Q9_rect.cpp"
#pragma once

/* Here, we define the density in the outlet (which is actually the *incoming* density, on the left-hand side), the density in the inlet (which is actually the *outgoing* density, on the right-hand
 * side), the pressure drop across the system, and the Reynolds number for the simulation. */
const double rho_outlet = 1.0;
double dpdx = 1.0e-15;
double rho_inlet = 1.0;
double Reynolds_num = (y_len * u_max_centre)/kinemat_visc_FECN;


/* Here, we define the array that has the analytic solution to the problem, and the array that has the error between the simulation and the analytic solution. For problems with analytic solutions, we
 * use this to calculate the L₂ error, which allows us to verify the accuracy of this simulation. */
vect1D_double ux_analytic(y_len);
vect1D_double error(x_len);


/* Here, we define the quantities that allow us to examine whether the x-velocities have converged. prev_ux_profile is used to store the x-velocity profile at the previous timestep, for comparison.
 * The simulation is considered to have converged when the largest difference in any *individual* ux (i.e., the largest difference between ux[i][j] at *any* point on the lattice) is less than the
 * convergence threshold. Here, this amounts to requiring largest_diff_in_simulated_ux_between_timesteps < convergence_threshold_for_simulated_ux_diff. */
double largest_diff_in_simulated_ux_between_timesteps = 100.0;
const double convergence_threshold_for_simulated_ux_diff = 1.0e-10;
array2D_double prev_ux_profile(x_len, vect1D_double(y_len, 0.0));


/* Here, we apply the macroscopic boundary conditions for 2D Poiseuille flow. As mentioned in Krüger Sec. 5.3.4.1 (Pg. 190), since the lattice Boltzmann equation doesn't directly deal with the
 * macroscopic fields (density and velocity), we need to explicitly impose the continuity equation, the no-slip condition, and the no-penetration condition.
 * 
 * We note that the macroscopic boundary conditions need to be defined *per the specific problem*, since we need to satisfy this per the *given* geometry and boundary conditions. Thus, we can't
 * define a *general* function to impose the macroscopic boundary conditions (unless we're doing something fancy like recreating Palabos); instead, we need to define the macroscopic boundary
 * condition function *per problem*. However, we have the advantage that this function gives us the macroscopic BCs for *all* rectangular channel problems, since this is purely in terms of the wall
 * velocities (given by u_wall_bottom and u_wall top).
 *
 * The wall density is straightforwardly calculated in the aforementioned section of Krüger. However, we note that there is an explicit de-linking of timestep and space-step in this derivation. Thus,
 * for a speed c = dx/dt for timestep dt and space step dx, expressions that should be c/(c ± u_wall) are instead written here as 1/(1 ± u_wall). This uses Eq. 5.32 and Eq. 5.33 (Pg. 191) in
 * Krüger (in Sec. 5.3.4.1). */
void macroscopic_BCs_Poiseuille(array2D_double &rho_simulated_on_lattice, array2D_double &ux_simulated_on_lattice, array2D_double &uy_simulated_on_lattice, array3D_double propagated_dist_fn,
                                                                          const double bottom_wall_vel, const double top_wall_vel) {
  int x_amt = rho_simulated_on_lattice.size();
  int y_amt = rho_simulated_on_lattice[0].size();

  for (int i = 0; i < x_amt; i++) {
    // Here, we apply the boundary conditions for the bottom wall: continuity, no-slip (ux_simulated_on_lattice[i][0] = u_wall_bottom), and no-penetration (uy_simulated_on_lattice[i][0] = 0.0).
    rho_simulated_on_lattice[i][0] = (1.0/(1.0 - uy_simulated_on_lattice[i][0])) *
           (propagated_dist_fn[i][0][0] + propagated_dist_fn[i][0][1] + propagated_dist_fn[i][0][3] + 2.0 * (propagated_dist_fn[i][0][4] + propagated_dist_fn[i][0][7] + propagated_dist_fn[i][0][8]));
    ux_simulated_on_lattice[i][0] = bottom_wall_vel;
    uy_simulated_on_lattice[i][0] = 0.0;

    // Here, we apply the boundary conditions for the top wall: continuity, no-slip (ux_simulated_on_lattice[i][0] = u_wall_top), and no-penetration (uy_simulated_on_lattice[i][0] = 0.0).
    rho_simulated_on_lattice[i][y_amt - 1] = (1.0/(1.0 + uy_simulated_on_lattice[i][y_amt - 1])) * (propagated_dist_fn[i][y_amt - 1][0] + propagated_dist_fn[i][y_amt - 1][1]
                                      + propagated_dist_fn[i][y_amt - 1][3] + 2.0 * (propagated_dist_fn[i][y_amt - 1][2] + propagated_dist_fn[i][y_amt - 1][5] + propagated_dist_fn[i][y_amt - 1][6]));
    ux_simulated_on_lattice[i][y_amt - 1] = top_wall_vel;
    uy_simulated_on_lattice[i][y_amt - 1] = 0.0;

  }
}


/* Here, we apply the pressure boundary conditions for 2D Poiseuille flow. As mentioned earlier, we have a halo region at x = 0 and x = x_len - 1; thus, we apply these boundary conditions at x = 1
 * and x = x_len - 2. We apply Eq. 7.30 in Krüger, on Pg. 286 (in Sec. 7.3.2). */
void pressure_BC_Poiseuille(array2D_double velocity_in_x_dir, array2D_double velocity_in_y_dir, vect1D_double x_dir_c_vals, vect1D_double y_dir_c_vals, vect1D_double lattice_dir_weights,
                                                              array3D_double current_equil_dist_fn, array3D_double &f_coll, array3D_double &u_dot_ci,
                                                              const double fluid_speed_of_sound_sq, const double inlet_density, const double outlet_density) {

  int x_dir_full = current_equil_dist_fn.size();
  int y_dir_full = current_equil_dist_fn[0].size();
  int q_directions = lattice_dir_weights.size();

  for (int j = 0; j < y_dir_full; j++) {
    for (int k = 0; k < q_directions; k++) {
      u_dot_ci[x_dir_full - 2][j][k] = velocity_in_x_dir[x_dir_full - 2][j] * x_dir_c_vals[k] + velocity_in_y_dir[x_dir_full - 2][j] * y_dir_c_vals[k];
      f_coll[0][j][k] = lattice_dir_weights[k] * (inlet_density + (u_dot_ci[x_dir_full - 2][j][k])/fluid_speed_of_sound_sq) + f_coll[x_dir_full - 2][j][k] - current_equil_dist_fn[x_dir_full - 2][j][k];

      u_dot_ci[1][j][k] = velocity_in_x_dir[1][j] * x_dir_c_vals[k] + velocity_in_y_dir[1][j] * y_dir_c_vals[k];
      f_coll[x_dir_full - 1][j][k] = lattice_dir_weights[k] * (outlet_density + u_dot_ci[1][j][k]/fluid_speed_of_sound_sq) + f_coll[1][j][k] - current_equil_dist_fn[1][j][k];
    }
  }

}


/* Here, we apply the equilibrium boundary conditions, which provide first-order accuracy, for 2D Couette flow. This uses the equilibrium scheme (ES), discussed in Sec. 5.3.4.2 (Pgs. 191-193) in
 * Krüger. Here, we use the expression given in Eq. 5.35 in Krüger (on Pg. 192, in Sec. 5.3.4.2). Since the pressure drop is zero here, that term goes away.
 * 
 * As with the macroscopic boundary conditions, we need to define the ES *per the specific problem*, since we need to satisfy this per the *given* geometry and boundary conditions. Thus, we can't
 * define a *general* function to impose the ES boundary conditions (unless we're doing something fancy like recreating Palabos); instead, we need to define the ES boundary condition function *per
 * problem*. */
void ES_Poiseuille(vect1D_double c_along_x_dir, vect1D_double c_along_y_dir, array2D_double x_dir_speed, array2D_double y_dir_speed, array2D_double dens_of_fluid, vect1D_double weight_of_lattice_dir,
                                                array3D_double &f_prop, array3D_double &u_dot_ci, const double speed_of_sound_sq) {
  int x_dir_qty = f_prop.size();
  int y_dir_qty = f_prop[0].size();
  int q_dir_qty = f_prop[0][0].size();

  for (int i = 0; i < x_dir_qty; i++) {
    for (int k = 0; k < q_dir_qty; k++) {
      u_dot_ci[i][0][k] = x_dir_speed[i][0] * c_along_x_dir[k] + y_dir_speed[i][0] * c_along_y_dir[k];
      u_dot_ci[i][y_dir_qty - 1][k] = x_dir_speed[i][y_dir_qty - 1] * c_along_x_dir[k] + y_dir_speed[i][y_dir_qty - 1] * c_along_y_dir[k];

      f_prop[i][0][k] = weight_of_lattice_dir[k] * (dens_of_fluid[i][0] + u_dot_ci[i][0][k]/speed_of_sound_sq);
      f_prop[i][y_dir_qty - 1][k] = weight_of_lattice_dir[k] * (dens_of_fluid[i][y_dir_qty - 1] + u_dot_ci[i][y_dir_qty - 1][k]/speed_of_sound_sq);
    
    }
  }
}


/* Here, we apply nonequilibrium boundary conditions, which provide second-order accuracy, for 2D Couette flow. This function in particular uses the Zou-He (nonequilibrium bounce back) boundary
 * conditions, discussed in Krüger Sec. 5.3.4.4, Pgs. 196-199. As with the equilibrium scheme and NEEM boundary conditions, we need to define the Zou-He boundary conditions *per the specific
 * problem*, since we need to satisfy this per the *given* geometry and boundary conditions. Thus, we can't define a *general* function to impose the Zou-He boundary conditions (unless we're doing
 * something fancy like recreating Palabos); instead, we need to define the Zou-He boundary condition function *per problem.* Note that, in particular, the Zou-He boundary conditions that are defined
 * involve solving a specific set of linear equations at the boundaries; this must usually be done by hand per problem. */
void Zou_He_Poiseuille(array2D_double x_dir_velocity, array2D_double y_dir_velocity, array2D_double dens_of_fluid, array3D_double &f_prop) {
  int x_dir_size = f_prop.size();
  int y_dir_size = f_prop[0].size();

  for (int i = 0; i < x_dir_size; i++) {
    /* Here, we incorporate the Zou-He (nonequilibrium bounce back) method, discussed in Krüger Sec. 5.3.4.4, Pgs. 196-199, Eq. 5.42-5.48. This first applies the ZH boundary conditions on the
     * bottom plate. */
    f_prop[i][0][2] = f_prop[i][0][4] + (2.0 * dens_of_fluid[i][0] * y_dir_velocity[i][0])/3.0;
    f_prop[i][0][5] = f_prop[i][0][7] + (dens_of_fluid[i][0] * y_dir_velocity[i][0])/6.0 - (f_prop[i][0][1] - f_prop[i][0][3])/2.0 + x_dir_velocity[i][0]/2.0;
    f_prop[i][0][6] = f_prop[i][0][8] + (dens_of_fluid[i][0] * y_dir_velocity[i][0])/6.0 + (f_prop[i][0][1] - f_prop[i][0][3])/2.0 - x_dir_velocity[i][0]/2.0;

    /* Here, we incorporate the Zou-He (nonequilibrium bounce back) method, discussed in Krüger Sec. 5.3.4.4, Pgs. 196-199, Eq. 5.42-5.48. This now applies the ZH boundary conditions on the top
     * plate. */
    f_prop[i][y_dir_size - 1][4] = f_prop[i][y_dir_size - 1][2] - (2.0 * dens_of_fluid[i][y_dir_size - 1] * y_dir_velocity[i][y_dir_size - 1])/3.0;
    f_prop[i][y_dir_size - 1][7] = f_prop[i][y_dir_size - 1][5] - (dens_of_fluid[i][y_dir_size - 1] * y_dir_velocity[i][y_dir_size - 1])/6.0
                                                                + (f_prop[i][y_dir_size - 1][1] - f_prop[i][y_dir_size - 1][3])/2.0 - x_dir_velocity[i][y_dir_size - 1]/2.0;
    f_prop[i][y_dir_size - 1][8] = f_prop[i][y_dir_size - 1][6] - (dens_of_fluid[i][y_dir_size - 1] * y_dir_velocity[i][y_dir_size - 1])/6.0
                                                                - (f_prop[i][y_dir_size - 1][1] - f_prop[i][y_dir_size - 1][3])/2.0 + x_dir_velocity[i][y_dir_size - 1]/2.0;
  
  }
}


/* Here, we define the function that compares the largest difference in the value of ux between timesteps. */
double previous_vs_current_ux_simulated(array2D_double simulated_ux_on_lattice_at_previous_iter, array2D_double simulated_ux_on_lattice_at_current_iter,
                                                                                                 double &largest_diff_in_simulated_ux_between_timesteps) {
  
  int x_dir_amt = simulated_ux_on_lattice_at_current_iter.size();
  int y_dir_amt = simulated_ux_on_lattice_at_current_iter[0].size();

  largest_diff_in_simulated_ux_between_timesteps = abs(simulated_ux_on_lattice_at_current_iter[0][1] - simulated_ux_on_lattice_at_previous_iter[0][1]);

  for (int i = 1; i < x_dir_amt; i++) {
    for (int j = 1; j < y_dir_amt - 1; j++) {
      float curr_place_diff = abs(simulated_ux_on_lattice_at_current_iter[i][j] - simulated_ux_on_lattice_at_previous_iter[i][j]);

      if (curr_place_diff > largest_diff_in_simulated_ux_between_timesteps) {
        largest_diff_in_simulated_ux_between_timesteps = curr_place_diff;
      }

      else {
        continue;
      }

    }
  }

  return largest_diff_in_simulated_ux_between_timesteps;
}


/* Here, we populate the array that has the analytic solution to the problem. For the Couette flow problem, this problem has infinite symmetry in the x-direction, and the velocity depends *only* on
 * the profile in the y-direction. In particular, the analytic solution should be ux = ux(y) = -(dp/dx) * y * (y - y_len)/ν (where y is a variable), uy = 0, and rho = 1 (since this is an
 * incompressible problem in the analytic case). */
vect1D_double calculate_analytic_sol_Poiseuille(vect1D_double &ux_analytic, const int y_height, const double u_max) {
  int poiseuille_channel_span = ux_analytic.size();
  float channel_height = y_height - 1.0;
  
  for (int i = 0; i < poiseuille_channel_span; i++) {
    ux_analytic[i] = (-4.0 * u_max)/(std::pow(channel_height, 2)) * i * (i - channel_height);
  }

  return ux_analytic;
}


/* Here, we calculate the L₂ error for 2D Poiseuille flow problems in rectangular coordinates. */
double calculate_L2_error_2D_Poiseuille_rectangular(vect1D_double &error, array2D_double x_dir_vel_simulated, vect1D_double analytic_x_dir_vel) {
  int sim_size = error.size() - 1;
  int channel_size = x_dir_vel_simulated[0].size();
  
  for (int i = 1; i < sim_size; i++) {
    double running_sum_ux_sim_vs_analytic = 0.0;
    double running_sum_ux_analytic = 0.0;

    for (int j = 0; j < channel_size; j++) {
      running_sum_ux_sim_vs_analytic += std::pow(x_dir_vel_simulated[i][j] - analytic_x_dir_vel[j], 2);
      running_sum_ux_analytic += std::pow(analytic_x_dir_vel[j], 2);
    }

    error[i] = std::sqrt(running_sum_ux_sim_vs_analytic/running_sum_ux_analytic);
  }

  double L2_error = std::accumulate(error.begin(), error.end(), 0.0)/sim_size;
  return L2_error;
}