#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h               = 0.67;
  double OmegaB          = 0.05;
  double OmegaCDM        = 0.267;
  double OmegaK          = 0.0;
  double Neff            = 3.046;
  double TCMB            = 2.7255;

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
  double kpivot_mpc      = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
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
  rec.solve(-13.0, 0.0, 10000, true, true, true, -7.0, true); // More points due to abrupt changes in dtau/dx and ddtau/ddx

  // Print freeze-out abundance of free electrons
  rec.print_freeze_out_abundance();

  // Print optical debths at reionization
  rec.print_tau_reionization(true);

  // Print the times and horizon sizes at decoupling and recombination (with and without Peebles)
  rec.print_decoupling_and_recombination(true);
  rec.print_decoupling_and_recombination(true, true);

  // Output recombination quantities
  rec.output(-12.0, 0.0, "results/recombination.txt", true, true, true, true);

  
  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // // Solve the perturbations
  // Perturbations pert(&cosmo, &rec);
  // pert.solve();
  // pert.info();
  
  // // Output perturbation quantities
  // double kvalue = 0.01 / Constants.Mpc;
  // pert.output(kvalue, "perturbations_k0.01.txt");
  
  // // Remove when module is completed
  // return 0;
  
  // //=========================================================================
  // // Module IV
  // //=========================================================================

  // PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  // power.solve();
  // power.output("cells.txt");
  
  // // Remove when module is completed
  // return 0;

  // Utils::EndTiming("Everything");
}
