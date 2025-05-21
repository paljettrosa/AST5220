#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  if (argc != 1 && argc != 4 && argc != 5) {
    // <fiducial>, <master> and <toy> are 1 (0) if (not) solving for fiducial, fiducial w/o Helium, reionization, polarization and neutrinos, and toy cosmology, respectively
    // <contrib> (optional) is 1 (0) if (not) computing the individual contributions to C_ell. Ignored if fiducial = 0
    // If no arguments are passed, only fiducial is solved for
    std::cout << "Usage: " << argv[0] << " <fiducial> <master> <toy> (<contrib>)" << std::endl;
    return 1;
  }

  int fiducial = 1;
  int master   = 0;
  int toy      = 0;
  int contrib  = 0;
  if (argc > 1) {
    fiducial = atoi(argv[1]);
    master   = atoi(argv[2]);
    toy      = atoi(argv[3]);
    if (argc == 5)
      contrib = atoi(argv[4]);
  }

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
  double Yp              = 0.245;
  double z_reion         = 8.0;
  double Delta_z_reion   = 0.5;
  double z_Hereion       = 3.5;
  double Delta_z_Hereion = 0.5;

  // Perturbation parameters
  int n_ell_Theta        = 11; 
  int n_ell_Theta_P      = 11;
  int n_ell_Nu           = 13;
  double k_min           = 0.00001/Constants.Mpc; 
  double k_max           = 1.0/Constants.Mpc;
  bool lensing           = true;

  // Power-spectrum parameters
  double A_s             = 2.1e-9; 
  double n_s             = 0.965; 
  double kpivot_Mpc      = 0.05;
  int ell_max            = 2500; 

  //=========================================================================
  // Fiducial cosmology
  //=========================================================================
  if (fiducial) {
    Utils::StartTiming("fiducial");

    //=========================================================================
    // Milestone I
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


    //=========================================================================
    // Milestone II
    //=========================================================================
    
    // Initialize instance of the recombination history class and print info
    RecombinationHistory rec(&cosmo, Yp, z_reion, Delta_z_reion, z_Hereion, Delta_z_Hereion);
    rec.info();

    // // Solve using only the Saha approximation
    // std::cout << "\nSolving for X_e and n_e using only the Saha approximation:\n";
    // rec.solve(-18.0, 0.0, 1000, 0.9999, false, false, false, -7.0, false, true); 
    
    // // Output recombination quantities
    // rec.output(-12.0, 0.0, "results/recombination_Saha.txt", false);

    // Solve using both Saha and Peebles
    std::cout << "\nSolving entire system with Saha and Peebles:\n";
    rec.solve(-18.0, 0.0, 10000, 0.9999, true, true, true, -7.0, true);

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
    // Milestone III
    //=========================================================================
  
    // Solve the perturbations
    Perturbations pert(&cosmo, &rec, n_ell_Theta, n_ell_Theta_P, n_ell_Nu, k_min, k_max, lensing);
    // Perturbations pert(&cosmo, &rec, n_ell_Theta, n_ell_Theta_P, n_ell_Nu, k_min, 0.1/Constants.Mpc, lensing); //TODO remove
    pert.info();
    pert.solve(-18.0, 0.0, 1000, 200); 
    // pert.solve(-18.0, 0.0, 500, 100); //TODO remove
    
    Vector k_values = {0.001/Constants.Mpc, 0.01/Constants.Mpc, 0.1/Constants.Mpc, 1.0/Constants.Mpc};
    std::vector<std::string> k_strings = {"0.001", "0.01", "0.1", "1.0"};
    for (int i; i < 4; i++) {
      // Print the time when tight-coupling ends
      pert.print_tight_coupling_time(k_values[i]);

      // Print the time of horizon entry
      pert.print_horizon_entry_time(k_values[i]);

      // Output perturbation quantities
      pert.output(-18.0, 0.0, k_values[i], "results/perturbations_k" + k_strings[i] + ".txt");
    }
    

    //=========================================================================
    // Milestone IV
    //=========================================================================

    PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_Mpc, ell_max);
    power.info();
    power.solve(-18.0, 0.0, 3000, k_min, 0.3/Constants.Mpc, false, true, true, true); //TODO add true at the end 

    // Print equality-scale
    power.print_equality_scale();

    // Output power spectrum quantities
    std::vector<int> ells = {15, 100, 350, 850, 1350, 1850}; 
    for (auto ell: ells)
      power.output_transfer_functions(k_min, 0.3/Constants.Mpc, ell, "results/transfer_functions_ell" + std::to_string(ell) + ".txt"); 
    power.output_C_of_theta("results/C_of_theta.txt", true);
    power.output_C_ells("results/C_ells.txt", false, true);
    power.output_P_k(0.000001, 5.0, "results/P_k.txt");
    power.output_xi(1.0, 500.0, "results/xi.txt");

    Utils::EndTiming("fiducial");

    // Compute individual contributions to the TT power spectrum 
    if (contrib) {
      Utils::StartTiming("contributions");

      // Sachs-Wolfe term
      pert.solve(-18.0, 0.0, 1000, 200, true, false, false, false); 
      PowerSpectrum power_SW(&cosmo, &rec, &pert, A_s, n_s, kpivot_Mpc, ell_max);
      power_SW.solve(-18.0, 0.0, 3000, k_min, 0.3/Constants.Mpc, true); 
      power_SW.output_C_ells("results/C_ell_SW.txt", true);

      // Integrated Sachs-Wolfe term
      pert.solve(-18.0, 0.0, 1000, 200, false, true, false, false); 
      PowerSpectrum power_ISW(&cosmo, &rec, &pert, A_s, n_s, kpivot_Mpc, ell_max);
      power_ISW.solve(-18.0, 0.0, 3000, k_min, 0.3/Constants.Mpc, true);
      power_ISW.output_C_ells("results/C_ell_ISW.txt", true);

      // Doppler term
      pert.solve(-18.0, 0.0, 1000, 200, false, false, true, false); 
      PowerSpectrum power_Doppler(&cosmo, &rec, &pert, A_s, n_s, kpivot_Mpc, ell_max);
      power_Doppler.solve(-18.0, 0.0, 3000, k_min, 0.3/Constants.Mpc, true);
      power_Doppler.output_C_ells("results/C_ell_Doppler.txt", true);

      // Polarization term
      pert.solve(-18.0, 0.0, 1000, 200, false, false, false, true); 
      PowerSpectrum power_pol(&cosmo, &rec, &pert, A_s, n_s, kpivot_Mpc, ell_max);
      power_pol.solve(-18.0, 0.0, 3000, k_min, 0.3/Constants.Mpc, true);
      power_pol.output_C_ells("results/C_ell_polarization.txt", true);

      Utils::EndTiming("contributions");
    }
  }

  //=========================================================================
  // Fiducial cosmology w/o Helium, reionization, polarization or neutrinos
  //=========================================================================
  if (master) {
    Utils::StartTiming("master");

    BackgroundCosmology cosmo_master(h, Omega_b0, Omega_CDM0, Omega_k0, 0.0, T_CMB0);
    cosmo_master.solve(-21.0, 6.0, 1000);
    RecombinationHistory rec_master(&cosmo_master, 0.0, 0.0, 0.0, 0.0, 0.0);
    rec_master.solve(-18.0, 0.0, 10000, 0.99, true);
    Perturbations pert_master(&cosmo_master, &rec_master, n_ell_Theta, 0, 0, k_min, k_max, false);
    pert_master.solve(-18.0, 0.0, 1000, 200); 
    PowerSpectrum power_master(&cosmo_master, &rec_master, &pert_master, A_s, n_s, kpivot_Mpc, ell_max);
    power_master.solve(-18.0, 0.0, 3000, k_min, 0.3/Constants.Mpc, true); 
    power_master.output_C_ells("results/C_ell_master.txt", true);

    Utils::EndTiming("master");
  }

  //=========================================================================
  // Toy cosmology
  //=========================================================================
  if (toy) {
    Utils::StartTiming("toy");

    // Milestone I
    BackgroundCosmology cosmo_toy(0.7, 0.05, 0.45, 0.0, 0.0, T_CMB0);
    cosmo_toy.solve(-21.0, 6.0, 1000);
    cosmo_toy.output(-20.0, 5.0, "results/toy/cosmology.txt");

    // Milestone II
    RecombinationHistory rec_toy(&cosmo_toy, 0.0, 0.0, 0.0, 0.0, 0.0);
    rec_toy.solve(-18.0, 0.0, 1000, 0.99);
    rec_toy.output(-12.0, 0.0, "results/toy/recombination.txt");

    // Milestone II with reionization
    RecombinationHistory rec_reion(&cosmo_toy, 0.0, 11.0, Delta_z_reion, 0.0, 0.0);
    rec_reion.solve(-18.0, 0.0, 1000, 0.99);
    rec_reion.output(-12.0, 0.0, "results/toy/recombination_reion.txt");

    // Milestone II with Helium and reionization
    RecombinationHistory rec_Hereion(&cosmo_toy, 0.24, 11.0, Delta_z_reion, z_Hereion, Delta_z_Hereion);
    rec_Hereion.solve(-18.0, 0.0, 1000, 0.99);
    rec_Hereion.output(-12.0, 0.0, "results/toy/recombination_Hereion.txt");

    // Milestone III
    Perturbations pert_toy(&cosmo_toy, &rec_toy, 8, 0, 0, 0.00005/Constants.Mpc, 0.3/Constants.Mpc, false);
    pert_toy.solve(-18.0, 0.0, 1000, 200);
    pert_toy.output(-18.0, 0.0, 0.001/Constants.Mpc, "results/toy/perturbations_k0.001.txt");
    pert_toy.output(-18.0, 0.0, 0.01/Constants.Mpc, "results/toy/perturbations_k0.01.txt");
    pert_toy.output(-18.0, 0.0, 0.1/Constants.Mpc, "results/toy/perturbations_k0.1.txt");
    
    // Milestone IV
    PowerSpectrum power_toy(&cosmo_toy, &rec_toy, &pert_toy, 1e-9, n_s, kpivot_Mpc, 2000);
    power_toy.solve(-18.0, 0.0, 3000, 0.00005/Constants.Mpc, 0.3/Constants.Mpc);
    power_toy.output_C_ells("results/toy/C_TT.txt");
    power_toy.output_P_k(0.00005, 0.3, "results/toy/P_k.txt");

    Utils::EndTiming("toy");
  }
}
