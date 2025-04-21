#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h               = 0.67;
  double Omega_b0        = 0.05;
  double Omega_CDM0      = 0.267;
  double Omega_k0        = 0.0;
  double N_eff           = 3.046;
  double T_CMB0          = 2.7255;

  // Recombination/reionization parameters
  // double Yp              = 0.0;
  // double z_reion         = 0.0;
  // double Delta_z_reion   = 0.0;
  // double z_Hereion       = 0.0;
  // double Delta_z_Hereion = 0.0;
  double Yp              = 0.245;
  double z_reion         = 8.0;
  double Delta_z_reion   = 0.5;
  double z_Hereion       = 3.5;
  double Delta_z_Hereion = 0.5;
  // double z_Hereion       = 0.0;
  // double Delta_z_Hereion = 0.0;

  // Power-spectrum parameters
  double A_s             = 2.1e-9;
  double n_s             = 0.965;
  double kpivot_Mpc      = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, Omega_b0, Omega_CDM0, Omega_k0, N_eff, T_CMB0);
  cosmo.info();
  cosmo.solve(-21.0, 6.0, 1000, true, true);

  // Print the current day values of the cosmic and conformal times
  cosmo.print_times();
  
  // Output background evolution quantities
  cosmo.output(-20.0, 5.0, "results/cosmology.txt", true, true, true, true);

  // // Do the supernova fits
  // Vector bestfit_params {0.0, 0.0, 0.0};
  // mcmc_fit_to_supernova_data("data/supernovadata.txt", "results/results_supernovafitting.txt", &bestfit_params);

  // // Solve background with best-fit parameters
  // BackgroundCosmology cosmo_bestfit(bestfit_params[0], OmegaB, bestfit_params[1] - OmegaB, bestfit_params[2], Neff, TCMB);
  // cosmo_bestfit.info();
  // cosmo_bestfit.solve();
  
  // // Output background evolution quantities
  // cosmo_bestfit.output(-20.0, 5.0, "results/cosmology_bestfit.txt", false, false, true, false);

  // Remove when module is completed
  // return 0;

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Initialize instance of the recombination history class and print info
  RecombinationHistory rec(&cosmo, Yp, z_reion, Delta_z_reion, z_Hereion, Delta_z_Hereion);
  rec.info();

  // // Solve using only the Saha approximation
  // std::cout << "\nSolving for X_e and n_e using only the Saha approximation:\n";
  // rec.solve(-13.0, 0.0, 1000, false, false, false, -7.0, false, true); 
  
  // // Output recombination quantities
  // rec.output(-12.0, 0.0, "results/recombination_Saha.txt", false);

  // Solve using both Saha and Peebles
  std::cout << "\nSolving entire system with Saha and Peebles:\n";
  rec.solve(-18.0, 0.0, 10000, true, true, true, -7.0, true); // More points due to abrupt changes in dtau/dx and ddtau/ddx

  // // Print freeze-out abundance of free electrons
  // rec.print_freeze_out_abundance();

  // // Print optical debths at reionization
  // rec.print_tau_reionization(true);

  // // Print the times and horizon sizes at decoupling and recombination (with and without Peebles)
  // rec.print_decoupling_and_recombination(true);
  // rec.print_decoupling_and_recombination(true, true);

  // Output recombination quantities
  rec.output(-12.0, 0.0, "results/recombination.txt", true, true, true, true);

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  // Perturbations pert(&cosmo, &rec);
  // pert.solve(-18.0, 0.0, 10000, 0.00005 / Constants.Mpc, 0.3 / Constants.Mpc, 200);
  BackgroundCosmology cosmo_toy(0.7, 0.05, 0.45, 0.0, 0.0, T_CMB0);
  cosmo_toy.solve(-21.0, 6.0, 1000);
  cosmo_toy.output(-20.0, 5.0, "results/cosmology_toy.txt");
  RecombinationHistory rec_toy(&cosmo_toy, 0.0, 0.0, 0.0, 0.0, 0.0); //TODO: add back x_of_tau and change Xe_Saha_limit back
  rec_toy.solve(-18.0, 0.0, 10000);
  rec_toy.output(-18.0, 0.0, "results/recombination_toy.txt");

  Perturbations pert(&cosmo_toy, &rec_toy, false, false, false);
  pert.info();
  pert.solve(-18.0, 0.0, 10000, 0.00005/Constants.Mpc, 0.3/Constants.Mpc, 200);
  
  // Output perturbation quantities
  double kvalue1 = 0.001 / Constants.Mpc;
  double kvalue2 = 0.01 / Constants.Mpc;
  double kvalue3 = 0.1 / Constants.Mpc;
  pert.output(-18.0, 0.0, kvalue1, "results/perturbations_k0.001.txt");
  pert.output(-18.0, 0.0, kvalue2, "results/perturbations_k0.01.txt");
  pert.output(-18.0, 0.0, kvalue3, "results/perturbations_k0.1.txt");
  
  // // Remove when module is completed
  // return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo_toy, &rec_toy, &pert, A_s, n_s, kpivot_Mpc);
  power.info();
  power.solve(-18.0, 0.0, 1000, 0.00005/Constants.Mpc, 0.3/Constants.Mpc, 200);
  power.output_P_k(0.00005, 0.3, "results/P_k.txt");
  power.output_C_ells("results/C_ells.txt");

  Utils::EndTiming("everything");
}
