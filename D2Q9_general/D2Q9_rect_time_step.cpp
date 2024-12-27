/* This is the *.cpp file that defines all of the time steppers used for the square-coordinate D2Q9 lattice Boltzmann (LB) simulations. All of these simulations are isothermal.
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
#include "./D2Q9_rect.cpp"
#pragma once


/* Here, we define the BGK collision operator, which we use in subsequent functions for the forward time stepping. */
array3D_double BGK_collision(array3D_double dist_fn_at_previous_time, array3D_double eq_dist_fn_at_previous_time, const double BGK_relaxation_time) {
  int x_dir_span = dist_fn_at_previous_time.size();
  int y_dir_span = dist_fn_at_previous_time[0].size();
  int q_dir_span = dist_fn_at_previous_time[0][0].size();
  
  array3D_double BGK_collision_term = eq_dist_fn_at_previous_time;

  for (int i = 0; i < x_dir_span; i++) {
    for (int j = 0; j < y_dir_span; j++) {
      for (int k = 0; k < q_dir_span; k++) {
        BGK_collision_term[i][j][k] = - (1.0/BGK_relaxation_time) * (dist_fn_at_previous_time[i][j][k] - eq_dist_fn_at_previous_time[i][j][k]);
      }
    }
  }

  return BGK_collision_term;
}
                                                                                        
/* Here, we perform the collision step, using the "standard" lattice Boltzmann update expression. This uses Eq. 3.13 in Krüger Sec. 3.3.3.5 (Pg. 70). Across the collision, streaming, and boundary
 * condition steps, the collision and propagated distribution functions (f_coll and f_prop respectively here, f_i and f_i^* in Krüger Ch. 4) swap places (during the collision step) and then swap
 * places again (during the streaming step.) Here, we collide (update f_coll). We will then stream (update f_prop) and then apply boundary conditions to f_prop. 
 * 
 * As mentioned in Krüger Sec. 3.5.2 (Pgs. 96-98), the "standard" lattice Boltzmann expression is a 1st-order forward Euler time stepper. However, the 2nd-order Crank-Nicolson (trapezoidal) time
 * stepper looks identical to the 1st-order forward Euler time stepper, since the two are related by a variable transformation. This is mentioned in that section, and explained very well in D. Wilde
 * et al. Thus, this time stepper is 2nd-order accurate. The BGK_collision function used here is defined on line 42 of this *.cpp file. */
array3D_double collision_forward_Eul_Crank_Nic(array3D_double &f_coll, array3D_double previously_propagated_dist_fn, array3D_double current_eq_dist_fn, const double BGK_relax_time) {
  int x_span = f_coll.size();
  int y_span = f_coll[0].size();
  int q_span = f_coll[0][0].size();

  array3D_double collision_term = BGK_collision(previously_propagated_dist_fn, current_eq_dist_fn, BGK_relax_time);

  for (int i = 0; i < x_span; i++) {
    for (int j = 0; j < y_span; j++) {
      for (int k = 0; k < q_span; k++) {
        f_coll[i][j][k] = previously_propagated_dist_fn[i][j][k] + collision_term[i][j][k];

        //We could also have used inverse_relaxation_time = 1.0/BGK_relaxation_time. Defining the input "omega" to stand in for inverse_relaxation_time, this would have appeared as:
        //f_coll[i][j][k] = (1 - omega)) * previously_propagated_dist_fn[i][j][k] + omega * current_eq_dist_fn[i][j][k];
      }
    }
  }

  return f_coll;
}


/* Here, we perform the collision step, using the 3rd-order accurate Adams-Moulton (AM3) timestepper provided in D. Wilde et al. The paper discusses how to create explicit multistep time steppers out
 * of implicit methods (such as AM3). Here, we use the explicit formula for AM3 time stepping given as Eq. 11 in D. Wilde et al. (on Pg. 157; the fourth page in the article). The BGK_collision
 * function used here is defined on line 42 of this *.cpp file. */
array3D_double collision_AM3_Wilde(array3D_double &f_coll, array3D_double previously_propagated_dist_fn_current_timestep, array3D_double current_timestep_eq_dist_fn, 
                                                           array3D_double previously_propagated_dist_fn_previous_timestep, array3D_double previous_timestep_eq_dist_fn,
                                                           const double BGK_relax_time_current_timestep, const double BGK_relax_time_previous_timestep) {
  int x_span = f_coll.size();
  int y_span = f_coll[0].size();
  int q_span = f_coll[0][0].size();

  array3D_double collision_current = BGK_collision(previously_propagated_dist_fn_current_timestep, current_timestep_eq_dist_fn, BGK_relax_time_current_timestep);
  array3D_double collision_previous = BGK_collision(previously_propagated_dist_fn_previous_timestep, previous_timestep_eq_dist_fn, BGK_relax_time_previous_timestep);

  for (int i = 0; i < x_span; i++) {
    for (int j = 0; j < y_span; j++) {
      for (int k = 0; k < q_span; k++) {
        f_coll[i][j][k] = previously_propagated_dist_fn_current_timestep[i][j][k] + collision_current[i][j][k] - collision_previous[i][j][k];

        //f_coll[i][j][k] = previously_propagated_dist_fn_current_timestep[i][j][k] 
        //                                            - (1.0/BGK_relax_time_current_timestep) * (previously_propagated_dist_fn_current_timestep[i][j][k] - current_eq_dist_fn[i][j][k])
        //                                            + (1.0/BGK_relax_time_previous_timestep) * (previously_propagated_dist_fn_previous_timestep[i][j][k] - previous_timestep_eq_dist_fn[i][j][k]);
      }
    }
  }

  return f_coll;
}


/* Here, we perform the collision step, using the 2nd-order accurate backwards differentiation formula (BDF2) timestepper provided in D. Wilde et al. The paper discusses how to create explicit
 * multistep time steppers out of implicit methods (such as BDF2). Here, we use the explicit formula for BDF2 time stepping given as Eq. 12 in D. Wilde et al. (on Pg. 157; the fourth page in the
 * article). The BGK_collision function used here is defined on line 42 of this *.cpp file. */
array3D_double collision_BDF2_Wilde(array3D_double &f_coll, array3D_double previously_propagated_dist_fn_current_timestep, array3D_double current_timestep_eq_dist_fn,
                                                            array3D_double previously_propagated_dist_fn_previous_timestep, array3D_double previous_timestep_eq_dist_fn,
                                                            const double BGK_relax_time_current_timestep, const double BGK_relax_time_previous_timestep) {
  int x_span = f_coll.size();
  int y_span = f_coll[0].size();
  int q_span = f_coll[0][0].size();

  array3D_double collision_current = BGK_collision(previously_propagated_dist_fn_current_timestep, current_timestep_eq_dist_fn, BGK_relax_time_current_timestep);
  array3D_double collision_previous = BGK_collision(previously_propagated_dist_fn_previous_timestep, previous_timestep_eq_dist_fn, BGK_relax_time_previous_timestep);

  for (int i = 0; i < x_span; i++) {
    for (int j = 0; j < y_span; j++) {
      for (int k = 0; k < q_span; k++) {
        f_coll[i][j][k] = (4.0/3.0) * previously_propagated_dist_fn_current_timestep[i][j][k] - previously_propagated_dist_fn_previous_timestep[i][j][k]/3.0
                                      + collision_current[i][j][k] - collision_previous[i][j][k];

        //f_coll[i][j][k] = (4.0/3.0) * previously_propagated_dist_fn_current_timestep[i][j][k] - previously_propagated_dist_fn_previous_timestep[i][j][k]/3.0
        //                                            - (1.0/BGK_relax_time_current_timestep) * (previously_propagated_dist_fn_current_timestep[i][j][k] - current_timestep_eq_dist_fn[i][j][k])
        //                                            + (1.0/BGK_relax_time_previous_timestep) * (previously_propagated_dist_fn_previous_timestep[i][j][k] - previous_timestep_eq_dist_fn[i][j][k]);
      }
    }
  }

  return f_coll;
}


/* Here, we perform the collision step, using a "standard" explicit Runge-Kutta 4th order time stepper. As mentioned in Krüger Sec. 3.5.2 (Pgs. 96-98), the "standard" lattice Boltzmann discretisation
 * involves coupling time and space together. As such, to use a scheme like RK4 or TVD-RK3, we generally need to decouple the time and space discretisations. This can be done using the planned finite
 * volume LBM (FV-LBM) approach, reproducing DUGKS, and/or an explicit discretisation of space (as in the Shirasat-Nayak-Patil paper). These are explicitly the planned extensions for this code, so
 * this isn't an impediment (and, indeed, this function needs to be checked when those extensions are performed). The BGK_collision function used here is defined on line 42 of this *.cpp file. */
array3D_double collision_RK4(array3D_double &f_coll, array3D_double previously_propagated_dist_fn, array3D_double current_eq_dist_fn, const double relaxation_time) {
  int x_span = f_coll.size();
  int y_span = f_coll[0].size();
  int q_span = f_coll[0][0].size();

  /* Here, we set the intermediate collision operators. The collision term for the predictor step we can calculate immediately. All of the other terms we can simply initialise right now; they will
   * have to be calculated after each step is performed. */
  array3D_double collision_term_intermediate_predict = BGK_collision(previously_propagated_dist_fn, current_eq_dist_fn, relaxation_time); // This is k_1 in the standard RK4 expression.

  array3D_double collision_term_intermediate_correct_1st = current_eq_dist_fn; // This will be k_2 in the standard RK4 expression.
  array3D_double collision_term_intermediate_correct_2nd = current_eq_dist_fn; // This will be k_3 in the standard RK4 expression.
  array3D_double collision_term_intermediate_correct_3rd = current_eq_dist_fn; // This will be k_4 in the standard RK4 expression.

  array3D_double f_prop_intermediate = current_eq_dist_fn; // This is the "intermediate stage" f_prop.

  /* Here, we perform the first (predictor) step.  */
  for (int i = 0; i < x_span; i++) {
    for (int j = 0; j < y_span; j++) {
      for (int k = 0; k < q_span; k++) {
        f_prop_intermediate[i][j][k] = previously_propagated_dist_fn[i][j][k] + collision_term_intermediate_predict[i][j][k];
      }
    }
  }
  collision_term_intermediate_correct_1st = BGK_collision(f_prop_intermediate, current_eq_dist_fn, relaxation_time);

  /* Here, we perform the first corrector step. */
  for (int i = 0; i < x_span; i++) {
    for (int j = 0; j < y_span; j++) {
      for (int k = 0; k < q_span; k++) {
        f_prop_intermediate[i][j][k] = previously_propagated_dist_fn[i][j][k] + collision_term_intermediate_correct_1st[i][j][k];
      }
    }
  }
  collision_term_intermediate_correct_2nd = BGK_collision(f_prop_intermediate, current_eq_dist_fn, relaxation_time);

  /* Here, we perform the third corrector step. */
  for (int i = 0; i < x_span; i++) {
    for (int j = 0; j < y_span; j++) {
      for (int k = 0; k < q_span; k++) {
        f_prop_intermediate[i][j][k] = previously_propagated_dist_fn[i][j][k] + collision_term_intermediate_correct_2nd[i][j][k];
      }
    }
  }
  collision_term_intermediate_correct_3rd = BGK_collision(f_prop_intermediate, current_eq_dist_fn, relaxation_time);

  /* Here, we perform the final RK4 estimation. */
  for (int i = 0; i < x_span; i++) {
    for (int j = 0; j < y_span; j++) {
      for (int k = 0; k < q_span; k++) {
        f_coll[i][j][k] = previously_propagated_dist_fn[i][j][k]
                                              + (1.0/6.0) * (collision_term_intermediate_predict[i][j][k] + 2.0 * collision_term_intermediate_correct_1st[i][j][k]
                                                                 + 2.0 * collision_term_intermediate_correct_2nd[i][j][k] + collision_term_intermediate_correct_3rd[i][j][k]);
      }
    }
  }

  return f_coll;
}
