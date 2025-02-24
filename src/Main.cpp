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
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // Recombination parameters
  double Yp          = 0.245;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.info();
  cosmo.solve(true, true);

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
  return 0;

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  
  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.01.txt");
  
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
