/* This is the *.cpp file that carries out a simulation of 2D rectangular Coeutte flow using the compressible D2Q9 lattice Boltzmann (LB) scheme, over 17 different points in a single decade (order of
 * magnitude). This includes D2Q9_rect.cpp (which contains typedefs and functions for all simulations (compressible and incompressible alike), except the time steppers); D2Q9_rect_time_step.cpp
 * (which contains the time steppers), D2Q9_rect_compress.cpp (which contains functions used for compressible LB simulations), and D2Q9_Couette.cpp (which contains functions used for 2D Couette flow
 * problems specifically). */


/* Here, we include all of the packages necessary for the code. We also include D2Q9_rect.cpp, D2Q9_rect_time_step.cpp, D2Q9_rect_compress.cpp, and D2Q9_Couette.cpp. */
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
#include "./../../D2Q9_general/D2Q9_rect.cpp"
#include "./../../D2Q9_general/D2Q9_rect_time_step.cpp"
#include "./../../D2Q9_general/D2Q9_rect_compress.cpp"
#include "./../../D2Q9_Couette_Poiseuille/D2Q9_Couette.cpp"


/* Here, we finally begin the simulation. */
int main() {

  /* Here, we begin calculating the time that the simulation takes. */
  auto program_start = std::chrono::high_resolution_clock::now();

  /* Here, we begin saving the Reynolds number vs. L₂ error as *.CSVs. */
  std::ofstream reynolds_num_vs_L2_error_csv("./couette_flow_compressible_RK4/Re_vs_L2__u_wall_E5_by_1s.csv");
  reynolds_num_vs_L2_error_csv << "Re,L2" << std::endl;
  reynolds_num_vs_L2_error_csv << std::endl;

  /* Here, we run over each point per decade (17, to go from u_wall_top = 1.0e(OOM) to u_wall_top = 19.5e(OOM); where OOM is the order of magnitude set by hand). */
  for (int current_pt_in_decade = 0; current_pt_in_decade < 17; current_pt_in_decade++) {

    /* Here, we update the value of the top wall velocity, and calculate the corresponding Reynolds number. */
    u_wall_top = (u_wall_top_reference + current_pt_in_decade/2.0) * 1.0e5;
    Reynolds_num = (y_len * u_wall_top)/kinemat_visc_RK4;

    /* Here, we begin saving the x-velocities and densities as *.CSVs. */
    std::string Reynolds_num_string_for_filename = std::to_string(Reynolds_num);
    std::string ux_csv_path = "./couette_flow_compressible_RK4/ux_at_Re/";
    std::string ux_csv_filename = "ux_Re_" + Reynolds_num_string_for_filename;
    std::string rho_csv_path = "./couette_flow_compressible_RK4/rho_at_Re/";
    std::string rho_csv_filename = "rho_Re_" + Reynolds_num_string_for_filename;

    /* Here, we begin saving the values of the largest difference at each timestep for convergence. */
    //std::ofstream largest_diff_txt("./couette_flow_compressible_RK4/diff.txt");
    //largest_diff_txt << "Re = " << Reynolds_num << std::endl;


    ////////// THE INITIALISATION STARTS HERE //////////

    /* Here, we reset the macroscopic parameters (the density, x-velocity field, and y-velocity field). This allows us to start each simulation from a "fresh" set of macroscopic values, just in case.
     * This uses the reset_macroscopic_parameters function defined on line 162 of D2Q9_rect.cpp. */
    reset_macroscopic_parameters(rho_simulated_on_lattice, ux_simulated_on_lattice, uy_simulated_on_lattice, vel_sq_simulated_on_lattice, rho_global_start, ux_init_global, uy_init_global);

    /* Here, we initialise the equilibrium distribution function f_eq, as well as the distribution functions used in the collision step (f_coll) and the streaming step (f_prop). The process of
     * initialising the equilibrium distribution function is the same as the process of updating it at each timestep, using properly initialised expressions for the density and x- and y-velocity
     * arrays. These were properly initialised up front (when they were declared, on lines 123-130 of D2Q9_rect.cpp), so we can use the same compute_f_eq_compressible function as we use to compute
     * f_eq at each timestep. The vel_sq array may or may not be initialised properly, but since compute_f_eq_compressible updates vel_sq up front, it doesn't matter. Since this is a compressible
     * simulation, this uses the compute_f_eq_compressible function defined on line 53 of D2Q9_rect_compress.cpp. */
    compute_f_eq_compressible(f_eq, rho_simulated_on_lattice, ux_simulated_on_lattice, uy_simulated_on_lattice, vel_sq_simulated_on_lattice, cs_sq);
    f_coll = f_eq;
    f_prop = f_eq;


    ////////// THE MAIN SIMULATION LOOP STARTS HERE //////////

    /* Here, we perform the main simulation loop. Here, we follow the structure of the MATLAB scripts for Couette flow (couette_wetnode.m) and Poiseuille flow (poiseuille_wetnode.m), written by
     * Goncalo Silva as part of the code set accompanying the Krüger textbook. This script is located at https://github.com/lbm-principles-practice/code/blob/master/chapter5/couette_wetnode.m. In
     * particular, we use an order of the functions used largely based on that script, given by:
     *
     *            calculating macroscopic parameters -> imposing macroscopic BCs -> computing f_eq and colliding -> streaming -> imposing 1st order BCs -> imposing noneq BCs -> repeat
     *
     * The order of the functions used in couette_wetnode.m and poiseuille_wetnode.m is slightly different, where they impose *either* 1st order *or* noneq BCs, not both. The idea of imposing both is
     * inspired by the MATLAB script for the lid-driven cavity problem (cavity.m), written by Adriano Sciacovelli and packaged in Jonas Lätt's compilation of LBM script examples. A discussion of that
     * script is located in the Palabos forum post https://palabos-forum.unige.ch/t/lid-cavity-flow/31, mirrored at
     * https://web.archive.org/web/20230402025710/https://palabos-forum.unige.ch/t/lid-cavity-flow/31. We note that in addition to imposing both 1st order and noneq BCs, that script provides a
     * different order of functions used:
     *
     *            calculating macroscopic parameters -> imposing macroscopic BCs -> imposing noneq BCs -> computing f_eq and colliding -> imposing 1st order BCs -> streaming -> repeat
     *
     * cavity.m is mirrored at https://web.archive.org/web/20080616222518/http://www.lbmethod.org/_media/numerics:cavity.m. The order of operations used in cavity.m is discussed in detail in the
     * Palabos forum post mentioned earlier. Specifically, I took note of Jonas Lätt's comment: "Be aware of the fact that Zou/He boundary nodes need to execute a normal BGK collision. I know that
     * many people get this point wrong, and I therefore emphasize this once more: *Zou/He boundary nodes are required to execute a normal BGK collision step!* [emphasis original] In your code, this
     * simply means that implementation of Zou/He should come before collision [...]".
     *
     * The order of operations used in cavity.m is NOT the order of operations used in couette_wetnode.m or poiseuille_wetnode.m. We note that both orders of operations give accurate results for the
     * Couette flow problem. However, the order of operations used in couette_wetnode.m and poiseuille_wetnode.m (preferred by Krüger and Silva, which is followed here) gives an accurate result for
     * the Poiseuille flow problem. Conversely, the order of operations used in cavity.m (preferred by Sciacovelli and Lätt, which is not done here) does *not* give an accurate result for the
     * Poiseuille flow problem, but *does* give an accurate result here and for the cavity problem. */

    for (int t = 0; t < t_steps; t++) {
    //while (largest_diff_in_simulated_ux_between_timesteps > 1.0e-15) {

      /* Here, we compute the macroscopic quantities (density and the velocity fields) at each time step. Since this is a compressible simulation, this uses the
       * compute_macroscopic_params_compressible function defined on line 21 of D2Q9_rect_compress.cpp. */
      compute_macroscopic_params_compressible(f_prop, rho_simulated_on_lattice, ux_simulated_on_lattice, uy_simulated_on_lattice);

      /* Here, we apply the macroscopic boundary conditions. This uses the macroscopic_BCs_Couette function defined on line 40 of D2Q9_Couette.cpp. */
      macroscopic_BCs_Couette(rho_simulated_on_lattice, ux_simulated_on_lattice, uy_simulated_on_lattice, f_prop, u_wall_bottom, u_wall_top);

      /* Here, we compute the equilibrium distribution at the given timestep. Since this is a compressible simulation, this uses the compute_f_eq_compressible function defined on line 53 of
       * D2Q9_rect_compress.cpp. */
      compute_f_eq_compressible(f_eq, rho_simulated_on_lattice, ux_simulated_on_lattice, uy_simulated_on_lattice, vel_sq_simulated_on_lattice, cs_sq);

      /* Here, we perform the collision step. The "classical" scheme is a forward Euler (1st-order) scheme, although the Crank-Nicolson (trapezoidal 2nd-order) scheme is equivalent to this by a
       * variable transformation. This is performed via the collision_forward_Eul_Crank_Nic function, defined on line 61 of D2Q9_rect_time_step.cpp. The AM3 and BDF2 schemes are implicit schemes
       * that are turned into explicit schemes, following the technique in D. Wilde et al. The collision_AM3_Wilde function (for the 3rd-order Adams-Moulton method) is defined on line 90 of
       * D2Q9_rect_time_step.cpp, and the collision_BDF2_Wilde function (for the 2nd-order backwards differentiation formula) is defined on line 119 of D2Q9_rect_time_step.cpp. The RK4 time stepper
       * is homemade, using the general principle of RK4 time stepping. The collision_RK4 function that performs this is defined on line 149 of D2Q9_rect_time_step.cpp. */
      //collision_forward_Eul_Crank_Nic(f_coll, f_prop, f_eq, BGK_relax_time_current_FECN);
      //collision_AM3_Wilde(f_coll, f_prop, f_eq, f_prop_prev, f_eq_prev, BGK_relax_time_current_AM3, BGK_relax_time_prev_AM3);
      //collision_BDF2_Wilde(f_coll, f_prop, f_eq, f_prop_prev, f_eq_prev, BGK_relax_time_current_BDF2, BGK_relax_time_prev_BDF2);
      collision_RK4(f_coll, f_prop, f_eq, BGK_relax_time_current_RK4);

      /* Here, we perform the streaming step, without applying the boundary condition at the wall, but with applying the periodic boundary conditions. The boundary condition at the wall shifts
       * everything over, and this will be applied momentarily. This uses the streaming_periodic_x_and_y function defined on line 136 of D2Q9_rect.cpp. */
      streaming_periodic_x_and_y(cx_int, cy_int, f_prop, f_coll);

      /* Here, we apply the equilibrium conditions on the top and the bottom walls. The equilibrium boundary conditions provide first-order accuracy. We can either use the equil_BB_Couette function
       * defined on line 65 of D2Q9_Couette.cpp or the ES_Couette function defined on line 86 of D2Q9_Couette.cpp. */
      //equil_BB_Couette(u_wall_bottom, u_wall_top, cs_sq, D2Q9_lattice_weights, cx_float, rho, f_prop);
      ES_Couette(cx_float, cy_float, ux_simulated_on_lattice, uy_simulated_on_lattice, cs_sq, rho_simulated_on_lattice, D2Q9_lattice_weights, f_eq, f_prop, u_dot_ci);

      /* Here, we apply the nonequilibrium boundary conditions. The equilibrium boundary conditions provide first-order accuracy, whereas the nonequilibrium boundary conditions provide second-order
       * accuracy. We can use either the nonequilibrium extrapolation method (given here by NEEM_Couette, defined on line 118 of D2Q9_Couette.cpp) or the Zou-He function (given here by
       * Zou_He_Couette, defined on line 136 of D2Q9_Couette.cpp). Both are presented here, and since Couette flow has a linear profile, they both give the same results (and the same results as the
       * first-order boundary conditions). My personal preference is for the Zou-He method, both because it's the most intuitive to me, and because it's the most accurate amongst the three boundary
       * condition types (equilibrium, NEEM, and ZH). */
      //NEEM_Couette(cx_float, cy_float, ux_simulated_on_lattice, uy_simulated_on_lattice, cs_sq, rho_simulated_on_lattice, D2Q9_lattice_weights, f_eq, f_prop, u_dot_ci);
      Zou_He_Couette(ux_simulated_on_lattice, uy_simulated_on_lattice, rho_simulated_on_lattice, f_prop);

      /* Here, we compare the largest difference in values between the previous and current x-velocities. */
      previous_vs_current_ux_simulated(prev_ux_profile, ux_simulated_on_lattice, largest_diff_in_simulated_ux_between_timesteps);
      
      /* Here, we set the filename variables used to save the values of the largest difference at each timestep for convergence. */
      std::string time_step_string_for_filename = std::to_string(t);
      std::string Reynolds_num_string_for_filename = std::to_string(Reynolds_num);
      std::string curr_ux_csv_path = "./couette_flow_compressible_RK4/ux_profiles_for_testing_difference/";
      std::string curr_ux_csv_filename = "curr_ux_Re_" + Reynolds_num_string_for_filename + "__time_step_" + time_step_string_for_filename;
      std::string prev_ux_csv_filename = "prev_ux_Re_" + Reynolds_num_string_for_filename + "__time_step_" + time_step_string_for_filename;

      /* Here, we save the x-velocities, densities, and information about the largest difference in each timestep. We also save the current x-velocity profile for this part for the next timestep.
       * This uses the save_2D_array_to_csv function defined on line 197 of D2Q9_rect.cpp. */
      //save_2D_array_to_csv(ux_simulated_on_lattice, curr_ux_csv_path, curr_ux_csv_filename);
      //save_2D_array_to_csv(prev_ux_profile, curr_ux_csv_path, prev_ux_csv_filename);
      //std::ofstream largest_diff_txt("./couette_flow_compressible_RK4/diff.txt");
      //largest_diff_txt << largest_diff_in_simulated_ux_between_timesteps << std::endl;
      prev_ux_profile = ux_simulated_on_lattice;

    }

    //largest_diff_txt << std::endl << std::endl;

    /* Here, we calculate the L₂ error. This uses the calculate_analytic_sol_couette function defined on line 196 of D2Q9_Couette.cpp and the calculate_L2_error_2D_couette_poiseuille_rectangular
     * function defined on line 210 of D2Q9_Couette.cpp. */
    calculate_analytic_sol_Couette(ux_analytic, u_wall_top);
    double couette_sim_L2_error = calculate_L2_error_2D_Couette_rectangular(error, ux_simulated_on_lattice, ux_analytic);

    /* Here, we output the velocities and the densities to screen, to test against the analytic solution. This uses the print_2D_array function defined on line 182 of D2Q9_rect.cpp. */
    print_2D_array(ux_simulated_on_lattice);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    print_2D_array(rho_simulated_on_lattice);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    /* Here, we print the Reynolds number, L₂ error, and largest difference between x-velocity profiles at two different timesteps, all to screen. */
    std::cout << "Re = " << Reynolds_num << std::endl;
    std::cout << "L2 error: " << couette_sim_L2_error << std::endl;
    std::cout << "Largest iterated difference: " << largest_diff_in_simulated_ux_between_timesteps << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    /* Here, we save the Reynolds number vs. L₂ error, the densities, and the x-velocity profiles all to their respective *.CSVs. This uses the save_2D_array_to_csv function defined on line 197 of
     * D2Q9_rect.cpp. */
    reynolds_num_vs_L2_error_csv << Reynolds_num << "," << couette_sim_L2_error << std::endl;
    save_2D_array_to_csv(ux_simulated_on_lattice, ux_csv_path, ux_csv_filename);
    save_2D_array_to_csv(rho_simulated_on_lattice, rho_csv_path, rho_csv_filename);
  }

  /* Here, we calculate the amount of time the simulation took. */
  auto program_end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(program_end - program_start);

  /* Here, we print the amount of time the simulation took to screen. */
  std::cout << "Time taken by function: " << duration.count() << " seconds" << std::endl;
}
